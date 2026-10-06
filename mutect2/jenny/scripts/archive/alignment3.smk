import os
import sys
import yaml

# Create log directories
os.makedirs('logs/cluster/aln',exist_ok=True)
os.makedirs('logs/cluster/metrics',exist_ok=True)


### FUNCTIONS ###

def merge_input(wildcards):
    inputs=[]
    for _run in FILES[wildcards.sample].values():
        inputs.append(f"bam_input/work/{wildcards.sample}/{wildcards.reference}/{_run['run']}/{_run['lane']}/{_run['index']}/2.rg_sort.bam")
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
    def p(x):
        return x if x.startswith('FASTQ') or x.startswith('/') else 'FASTQ/'+x
    return {'R1':p(r1),'R2':p(r2)}

### ### PYTHON ### ###

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
            continue
        new_val={}
        for base_stem,r1,r2 in paired:
            parts=base_stem.split('_')
            run,lane,index=(parts[-3],parts[-2],parts[-1]) if len(parts)>=3 else ('run','0','0')
            new_val[base_stem]={'run':run,'lane':lane,'index':index,'files':[r1,r2]}
        FILES[sample]=new_val

with open(config['project']['sample_list']) as file:
    SAMPLES=file.read().splitlines()
for sample in SAMPLES:
    if sample not in FILES:
        print(f"WARNING: sample {sample} in sample.list not in fastq config, skipping",file=sys.stderr)
SAMPLES=[s for s in SAMPLES if s in FILES]
for sample in SAMPLES:
    os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

### ### ### RULES ### ### ###

localrules:aln_all,bam_table,write_bam_table

rule aln_all:
    input:
        expand("bam_input/final/{sample}/{sample}.{reference}.bam",sample=SAMPLES,reference=config['reference']['key']),
        expand("metrics/{reference}/{sample}/alignment_summary_metrics.txt",sample=SAMPLES,reference=config['reference']['key'])

rule bam_table:
    input:
        "bam.table"

rule bwa_mem:
    input:
        unpack(get_fastqs)
    output:
        temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"),
    log:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/bwa_mem.log"
    params:
        fasta=config['reference']['fasta'],
        outdir="bam_input/work/{sample}/{reference}/{run}/{lane}/{index}"
    threads:
        4
    resources:
        mem_mb=49152
    shell:
        """
        mkdir -p {params.outdir}
        (bwa-mem2 mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools view -bS -o {output}) &> {log}
        """

rule add_read_groups:
    input:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"
    output:
        temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/2.rg_sort.bam")
    log:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/add_read_groups.log"
    params:
        lb=config['resources']['library_key'],
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
            # --RGID {wildcards.sample}.{wildcards.run}.{wildcards.lane}.{wildcards.index}
            --RGPL illumina \
            --RGPU {wildcards.run}.{wildcards.lane}.{wildcards.index} \
            --RGSM {wildcards.sample} \
            --SORT_ORDER coordinate
        ) &> {log}
        """

rule merge_lanes:
    input:
        merge_input
    output:
        temp("bam_input/work/{sample}/{reference}/merged.bam")
    log:
        "bam_input/work/{sample}/{reference}/merge_lanes.log"
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
        "bam_input/work/{sample}/{reference}/merged.bam"
    output:
        temp("bam_input/work/{sample}/{reference}/input.bam")
    log:
        "bam_input/work/{sample}/{reference}/markdup.log"
    params:
        stats="bam_input/work/{sample}/{reference}/markdup.stats.txt",
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        (java -Xmx{params.memory} -jar $PICARDJAR MarkDuplicates \
            --INPUT {input} \
            --OUTPUT {output} \
            --METRICS_FILE {params.stats}
        ) &> {log}
        """


rule ValidateSamFile:
    input:
        "bam_input/work/{sample}/{reference}/input.bam"
    output:
        "metrics/{reference}/{sample}/validation_data.table"
    log:
        "bam_input/work/{sample}/{reference}/ValidateSamFile.log"
    params:
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        set +e
        mkdir -p metrics/{wildcards.reference}/{wildcards.sample}

        java -Xmx{params.memory} -jar $PICARDJAR ValidateSamFile I={input} O={output} MODE=SUMMARY 2> {log}
        exitcode=$?
        if [ $exitcode -ne 0 ]; then exit 1; fi
        exit 0
        """

rule CollectAlignmentSummaryMetrics:
    # Actual alignment QC (% mapped, mismatch rate, proper-pair rate) -- distinct from
    # ValidateSamFile above, which only checks that the BAM is well-formed, not how well
    # the reads aligned.
    input:
        "bam_input/work/{sample}/{reference}/input.bam"
    output:
        "metrics/{reference}/{sample}/alignment_summary_metrics.txt"
    log:
        "bam_input/work/{sample}/{reference}/CollectAlignmentSummaryMetrics.log"
    params:
        fasta=config['reference']['fasta'],
        memory="10240m"
    resources:
        mem_mb=14336
    shell:
        """
        mkdir -p metrics/{wildcards.reference}/{wildcards.sample}
        gatk --java-options -Xmx{params.memory} CollectAlignmentSummaryMetrics -R {params.fasta} -I {input} -O {output} &> {log}
        """

rule ready_bam:
    input:
        bam="bam_input/work/{sample}/{reference}/input.bam",
        table="metrics/{reference}/{sample}/validation_data.table"
    output:
        "bam_input/final/{sample}/{sample}.{reference}.bam"
    log:
        "bam_input/final/{sample}/{sample}.{reference}.ready_bam.log"
    shell:
        """
        (rsync -v {input.bam} {output}
        samtools index {output}) &> {log}
        """

rule write_bam_table:
    input:
        expand("bam_input/final/{sample}/{sample}.{reference}.bam",sample=SAMPLES,reference=config['reference']['key'])
    output:
        "input/bam.table"
    params:
        ref=config['reference']['key']
    run:
        with open(output[0],'w') as file:
            for sample in SAMPLES:
                file.write(f"{sample}\tbam_input/final/{sample}/{sample}.{params.ref}.bam\n")
