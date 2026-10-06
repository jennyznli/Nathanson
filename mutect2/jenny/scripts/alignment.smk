import os
import sys
import yaml

# Create log directories
os.makedirs('logs/cluster/aln',exist_ok=True)
os.makedirs('logs/cluster/metrics',exist_ok=True)

wildcard_constraints:
    sample="[^_/.]+",
    lb="[^/.]+",
    run="[^/.]+",
    lane="[^/.]+",
    index="[^/.]+",
    fq="[^/]+"

### FUNCTIONS ###

def p(x):
    return x if x.startswith('FASTQ') or x.startswith('/') else 'FASTQ/'+x

def bam_stem(sample):
    return f"alignment/{sample}/{sample}_{LIBRARY[sample]}"

def merge_input(wildcards):
    inputs=[]
    for _run in FILES[wildcards.sample].values():
        inputs.append(f"alignment/{wildcards.sample}/{wildcards.sample}.{_run['run']}.{_run['lane']}.{_run['index']}.rg.bam")
    assert len(inputs)>0
    return sorted(inputs)

def get_fastqs(wildcards):
    entry=None
    for _run in FILES[wildcards.sample].values():
        if _run['run']==wildcards.run and _run['lane']==wildcards.lane and _run['index']==wildcards.index:
            entry=_run
            break
    assert entry is not None
    r1,r2=entry['files'][0],entry['files'][1]
    return {'R1':p(r1),'R2':p(r2)}

### ### PYTHON ### ###

LIBRARY={}  # sample -> library field from FASTQ names (2nd '_' field), used for final BAM name + RGLB
FASTQS={}   # sample -> {fastq stem: path}, for FastQC

with open(config['project']['fastq_config']) as file:
    FILES=yaml.load(file,Loader=yaml.BaseLoader)
    for sample in list(FILES.keys()):
        val=FILES[sample]
        flist=sorted(val['files'])
        r1_by_base={}
        r2_by_base={}
        for f in flist:
            base=os.path.basename(f)
            if '_R1' in base:
                base_stem=base.replace('_R1.fastq.gz','').replace('_R1.fq.gz','')
                r1_by_base[base_stem]=f
            elif '_R2' in base:
                base_stem=base.replace('_R2.fastq.gz','').replace('_R2.fq.gz','')
                r2_by_base[base_stem]=f
        paired=[]
        for base_stem in sorted(set(r1_by_base)|set(r2_by_base)):
            r1=r1_by_base.get(base_stem)
            r2=r2_by_base.get(base_stem)
            if r1 is None:
                print(f"WARNING: unpaired R2 (no R1): {r2}",file=sys.stderr)
                continue
            if r2 is None:
                print(f"WARNING: unpaired R1 (no R2): {r1}",file=sys.stderr)
                continue
            paired.append((base_stem,r1,r2))
        if not paired:
            print(f"WARNING: no paired FASTQs for sample {sample}, skipping",file=sys.stderr)
            del FILES[sample]
            continue
        new_val={}
        libs=set()
        fqs={}
        for base_stem,r1,r2 in paired:
            parts=base_stem.split('_')
            run,lane,index=(parts[-3],parts[-2],parts[-1]) if len(parts)>=3 else ('run','0','0')
            new_val[base_stem]={'run':run,'lane':lane,'index':index,'files':[r1,r2]}
            libs.add(parts[1] if len(parts)>=2 else 'lib')
            for f in (r1,r2):
                fqs[os.path.basename(f).replace('.fastq.gz','').replace('.fq.gz','')]=p(f)
        # Final BAM is named {sample}_{library}, so a sample split across two library names is ambiguous
        if len(libs)!=1:
            sys.exit(f"ERROR: sample {sample} has FASTQs from multiple libraries ({','.join(sorted(libs))}); resolve before aligning")
        FILES[sample]=new_val
        LIBRARY[sample]=libs.pop()
        FASTQS[sample]=fqs

with open(config['project']['sample_list']) as file:
    SAMPLES=file.read().splitlines()
for sample in SAMPLES:
    if sample not in FILES:
        print(f"WARNING: sample {sample} in sample.list not in fastq config, skipping",file=sys.stderr)
SAMPLES=[s for s in SAMPLES if s in FILES]
for sample in SAMPLES:
    os.makedirs(f'logs/cluster/{sample}',exist_ok=True)
    os.makedirs(f"alignment/{sample}/logs",exist_ok=True)
    os.makedirs(f"alignment/{sample}/fastqc",exist_ok=True)
    os.makedirs(f"alignment/{sample}/mosdepth",exist_ok=True)

if not config['resources'].get('bqsr_known_sites'):
    sys.exit("ERROR: resources.bqsr_known_sites is empty in the config; BaseRecalibrator needs at least one known-sites VCF")

# Resource files are params (not inputs), so a dry run wouldn't notice them missing -- check up front
_resources=[config['reference']['fasta'],
            config['resources']['qc_targets_bed'],
            config['resources']['target_interval_list'],
            config['resources']['bait_interval_list'],
            config['resources']['bqsr_intervals']]
for vcf in config['resources']['bqsr_known_sites']:
    _resources+=[vcf,vcf+'.tbi']
_missing=[f for f in _resources if not os.path.exists(f)]
if _missing:
    sys.exit("ERROR: missing resource files (run process_targets_bed.sh / index known-sites VCFs?):\n  "+"\n  ".join(_missing))

QC_SUFFIXES=[
    'validation.txt',
    'markdup_metrics.txt',
    'recal_data.table',
    'alignment_summary_metrics.txt',
    'insert_size_metrics.txt',
    'hs_metrics.txt'
]

# Every per-sample QC file -- requested by aln_all and gathered by multiqc
ALL_QC=(
    [f"{bam_stem(s)}.{suffix}" for s in SAMPLES for suffix in QC_SUFFIXES]
    +[f"alignment/{s}/mosdepth/{s}_{LIBRARY[s]}.{m}.summary.txt" for s in SAMPLES for m in ('mosdepth','raw.mosdepth')]
    +[f"alignment/{s}/fastqc/{fq}_fastqc.zip" for s in SAMPLES for fq in FASTQS[s]]
)

### ### ### RULES ### ### ###

localrules:aln_all,bam_table,write_bam_table

rule aln_all:
    input:
        [f"{bam_stem(s)}.bam" for s in SAMPLES],
        ALL_QC,
        config['project']['bam_table'],
        "multiqc/multiqc_report.html"

rule bam_table:
    input:
        config['project']['bam_table']

rule fastqc:
    input:
        lambda wildcards: FASTQS[wildcards.sample][wildcards.fq]
    output:
        html="alignment/{sample}/fastqc/{fq}_fastqc.html",
        zip="alignment/{sample}/fastqc/{fq}_fastqc.zip"
    log:
        "alignment/{sample}/logs/{fq}.fastqc.log"
    threads:
        2
    shell:
        """
        fastqc --threads {threads} --outdir alignment/{wildcards.sample}/fastqc {input} &> {log}
        """

rule bwa_mem:
    input:
        unpack(get_fastqs)
    output:
        temp("alignment/{sample}/{sample}.{run}.{lane}.{index}.mapped.bam")
    log:
        "alignment/{sample}/logs/{sample}.{run}.{lane}.{index}.bwa_mem.log"
    params:
        fasta=config['reference']['fasta']
    threads:
        4
    resources:
        mem_mb=49152
    shell:
        """
        (bwa-mem2 mem -M -t {threads} \
            {params.fasta} \
            {input.R1} {input.R2} |
        samtools view -b -o {output}) &> {log}
        """

rule add_read_groups:
    input:
        "alignment/{sample}/{sample}.{run}.{lane}.{index}.mapped.bam"
    output:
        temp("alignment/{sample}/{sample}.{run}.{lane}.{index}.rg.bam")
    log:
        "alignment/{sample}/logs/{sample}.{run}.{lane}.{index}.add_read_groups.log"
    params:
        lb=lambda wildcards: f"{wildcards.sample}_{LIBRARY[wildcards.sample]}",
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        (java -Xmx{params.memory} -jar $PICARDJAR AddOrReplaceReadGroups \
            --INPUT {input} \
            --OUTPUT {output} \
            --RGID {wildcards.run}.{wildcards.lane} \
            --RGLB {params.lb} \
            --RGPL ILLUMINA \
            --RGPU {wildcards.run}.{wildcards.lane}.{wildcards.index} \
            --RGSM {wildcards.sample} \
            --SORT_ORDER coordinate
        ) &> {log}
        """

rule merge_lanes:
    input:
        merge_input
    output:
        temp("alignment/{sample}/{sample}.merged.bam")
    log:
        "alignment/{sample}/logs/{sample}.merge_lanes.log"
    params:
        inputs=lambda wildcards,input: " ".join(f"--INPUT {f}" for f in input),
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        (java -Xmx{params.memory} -jar $PICARDJAR MergeSamFiles \
            {params.inputs} \
            --OUTPUT {output} \
            --SORT_ORDER coordinate \
            --ASSUME_SORTED true \
            --USE_THREADING true
        ) &> {log}
        """

rule markdup:
    input:
        "alignment/{sample}/{sample}.merged.bam"
    output:
        bam=temp("alignment/{sample}/{sample}_{lb}.markdup.bam"),
        bai=temp("alignment/{sample}/{sample}_{lb}.markdup.bam.bai"),
        metrics="alignment/{sample}/{sample}_{lb}.markdup_metrics.txt"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.markdup.log"
    params:
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        (
        java -Xmx{params.memory} -jar $PICARDJAR MarkDuplicates \
            --INPUT {input} \
            --OUTPUT {output.bam} \
            --METRICS_FILE {output.metrics}

        samtools index {output.bam} {output.bai}
        ) &> {log}
        """

rule base_recalibrator:
    input:
        bam="alignment/{sample}/{sample}_{lb}.markdup.bam",
        bai="alignment/{sample}/{sample}_{lb}.markdup.bam.bai"
    output:
        "alignment/{sample}/{sample}_{lb}.recal_data.table"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.BaseRecalibrator.log"
    params:
        fasta=config['reference']['fasta'],
        intervals=config['resources']['bqsr_intervals'],
        known_sites=" ".join(f"--known-sites {f}" for f in config['resources']['bqsr_known_sites']),
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        gatk --java-options -Xmx{params.memory} BaseRecalibrator \
            -I {input.bam} \
            -R {params.fasta} \
            -L {params.intervals} \
            {params.known_sites} \
            -O {output} &> {log}
        """

rule apply_bqsr:
    input:
        bam="alignment/{sample}/{sample}_{lb}.markdup.bam",
        table="alignment/{sample}/{sample}_{lb}.recal_data.table"
    output:
        temp("alignment/{sample}/{sample}_{lb}.bqsr.bam")
    log:
        "alignment/{sample}/logs/{sample}_{lb}.ApplyBQSR.log"
    params:
        fasta=config['reference']['fasta'],
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        gatk --java-options -Xmx{params.memory} ApplyBQSR \
            -I {input.bam} \
            -R {params.fasta} \
            --bqsr-recal-file {input.table} \
            --create-output-bam-index false \
            -O {output} &> {log}
        """

rule validate_sam_file:
    input:
        "alignment/{sample}/{sample}_{lb}.bqsr.bam"
    output:
        "alignment/{sample}/{sample}_{lb}.validation.txt"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.ValidateSamFile.log"
    params:
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        set +e
        java -Xmx{params.memory} -jar $PICARDJAR ValidateSamFile I={input} O={output} MODE=SUMMARY 2> {log}
        exitcode=$?
        if [ $exitcode -ne 0 ]; then exit 1; fi
        exit 0
        """

rule ready_bam:
    input:
        bam="alignment/{sample}/{sample}_{lb}.bqsr.bam",
        table="alignment/{sample}/{sample}_{lb}.validation.txt"
    output:
        bam="alignment/{sample}/{sample}_{lb}.bam",
        bai="alignment/{sample}/{sample}_{lb}.bam.bai"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.ready_bam.log"
    shell:
        """
        (rsync -v {input.bam} {output.bam}
        samtools index {output.bam}) &> {log}
        """

### QC on the final BAM ###

rule collect_alignment_summary_metrics:
    input:
        "alignment/{sample}/{sample}_{lb}.bam"
    output:
        "alignment/{sample}/{sample}_{lb}.alignment_summary_metrics.txt"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.CollectAlignmentSummaryMetrics.log"
    params:
        fasta=config['reference']['fasta'],
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        gatk --java-options -Xmx{params.memory} CollectAlignmentSummaryMetrics \
            -I {input} \
            -O {output} \
            -R {params.fasta} &> {log}
        """

rule collect_insert_size_metrics:
    input:
        "alignment/{sample}/{sample}_{lb}.bam"
    output:
        metrics="alignment/{sample}/{sample}_{lb}.insert_size_metrics.txt",
        histogram="alignment/{sample}/{sample}_{lb}.insert_size_histogram.pdf"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.CollectInsertSizeMetrics.log"
    params:
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        gatk --java-options -Xmx{params.memory} CollectInsertSizeMetrics \
            -I {input} \
            -O {output.metrics} \
            -H {output.histogram} \
            -M 0.5 &> {log}
        """

rule collect_hs_metrics:
    input:
        "alignment/{sample}/{sample}_{lb}.bam"
    output:
        metrics="alignment/{sample}/{sample}_{lb}.hs_metrics.txt",
        per_target="alignment/{sample}/{sample}_{lb}.per_target_coverage.txt"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.CollectHsMetrics.log"
    params:
        fasta=config['reference']['fasta'],
        targets=config['resources']['target_interval_list'],
        baits=config['resources']['bait_interval_list'],
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        gatk --java-options -Xmx{params.memory} CollectHsMetrics \
            -I {input} \
            -O {output.metrics} \
            -R {params.fasta} \
            --BAIT_INTERVALS {params.baits} \
            --TARGET_INTERVALS {params.targets} \
            --PER_TARGET_COVERAGE {output.per_target} &> {log}
        """

rule mosdepth:
    input:
        bam="alignment/{sample}/{sample}_{lb}.bam",
        bai="alignment/{sample}/{sample}_{lb}.bam.bai"
    output:
        summary="alignment/{sample}/mosdepth/{sample}_{lb}.mosdepth.summary.txt",
        global_dist="alignment/{sample}/mosdepth/{sample}_{lb}.mosdepth.global.dist.txt",
        region_dist="alignment/{sample}/mosdepth/{sample}_{lb}.mosdepth.region.dist.txt",
        regions="alignment/{sample}/mosdepth/{sample}_{lb}.regions.bed.gz",
        thresholds="alignment/{sample}/mosdepth/{sample}_{lb}.thresholds.bed.gz",
        per_base="alignment/{sample}/mosdepth/{sample}_{lb}.per-base.bed.gz"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.mosdepth.log"
    params:
        prefix="alignment/{sample}/mosdepth/{sample}_{lb}",
        targets=config['resources']['qc_targets_bed']
    threads:
        4
    shell:
        """
        mosdepth \
            --threads {threads} \
            --by {params.targets} \
            --thresholds 50,100,250,500,1000 \
            {params.prefix} \
            {input.bam} &> {log}
        """

rule mosdepth_raw:
    input:
        bam="alignment/{sample}/{sample}_{lb}.bam",
        bai="alignment/{sample}/{sample}_{lb}.bam.bai"
    output:
        summary="alignment/{sample}/mosdepth/{sample}_{lb}.raw.mosdepth.summary.txt",
        global_dist="alignment/{sample}/mosdepth/{sample}_{lb}.raw.mosdepth.global.dist.txt",
        region_dist="alignment/{sample}/mosdepth/{sample}_{lb}.raw.mosdepth.region.dist.txt",
        regions="alignment/{sample}/mosdepth/{sample}_{lb}.raw.regions.bed.gz",
        thresholds="alignment/{sample}/mosdepth/{sample}_{lb}.raw.thresholds.bed.gz"
    log:
        "alignment/{sample}/logs/{sample}_{lb}.mosdepth_raw.log"
    params:
        prefix="alignment/{sample}/mosdepth/{sample}_{lb}.raw",
        targets=config['resources']['qc_targets_bed']
    threads:
        4
    shell:
        """
        mosdepth \
            --threads {threads} \
            --flag 772 \
            -n \
            --by {params.targets} \
            --thresholds 50,100,250,500,1000 \
            {params.prefix} \
            {input.bam} &> {log}
        """

rule write_bam_table:
    input:
        [f"{bam_stem(s)}.bam" for s in SAMPLES]
    output:
        config['project']['bam_table']
    run:
        with open(output[0],'w') as file:
            for sample in SAMPLES:
                file.write(f"{sample}\t{bam_stem(sample)}.bam\n")

rule multiqc:
    input:
        ALL_QC
    output:
        "multiqc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    resources:
        mem_mb=4000
    shell:
        """
        multiqc alignment --outdir multiqc --force &> {log}
        """
