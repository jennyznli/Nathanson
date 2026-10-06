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
        inputs.append(f"bam_input/work/{wildcards.sample}/{wildcards.reference}/{_run['run']}/{_run['lane']}/{_run['index']}/4.qsort.bam")
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
        # bwa-mem2's FM-index load alone wants ~28GB RAM (vs ~5GB for classic bwa) --
        # bumped from 32768 to leave headroom.
        mem_mb=49152
    shell:
        """
        mkdir -p {params.outdir}
        (bwa-mem2 mem -M -t {threads} {params.fasta} {input.R1} {input.R2} | samtools view -bS -o {output}) &> {log}
        """

rule samtools_readgroup:
    # Read-group tagging happens here, per lane, via addreplacerg -- not via bwa's -R --
    # so ID/PU (which differ per lane) and SM/LB (shared across a sample's lanes) come from
    # the actual wildcards instead of a hand-written, easy-to-typo/mis-escape -R string.
    input:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/1.mapped.bam"
    output:
        rg=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/2.readgroup.bam"),
        fm=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/3.fixmate.bam"),
        qs=temp("bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/4.qsort.bam"),
    log:
        "bam_input/work/{sample}/{reference}/{run}/{lane}/{index}/samtools_readgroup.log"
    params:
        LB=config['resources']['library_key']
    threads:
        4
    shell:
        """
        (
        samtools addreplacerg -@ {threads} -r 'ID:{wildcards.run}.{wildcards.lane}' -r 'PU:{wildcards.run}.{wildcards.lane}.{wildcards.index}' -r 'PL:illumina' -r 'LB:{params.LB}' -r 'SM:{wildcards.sample}' -o {output.rg} {input}
        samtools fixmate -m -@ {threads} {output.rg} {output.fm}
        samtools sort -@ {threads} -o {output.qs} {output.fm}
        ) &> {log}
        """

rule merge_lanes:
    # Merge every lane belonging to one sample BEFORE marking duplicates -- markdup can only
    # see duplicates within the file it's given, so per-lane dedup would miss PCR duplicates
    # of the same library that happened to land on a different sequencing lane.
    input:
        merge_input
    output:
        temp("bam_input/work/{sample}/{reference}/merged.bam")
    log:
        "bam_input/work/{sample}/{reference}/merge_lanes.log"
    threads:
        4
    shell:
        """
        samtools merge -f -@ {threads} {output} {input} &> {log}
        """

rule samtools_markdup:
    input:
        "bam_input/work/{sample}/{reference}/merged.bam"
    output:
        temp("bam_input/work/{sample}/{reference}/input.bam")
    log:
        "bam_input/work/{sample}/{reference}/samtools_markdup.log"
    params:
        stats="bam_input/work/{sample}/{reference}/markdup.stats.txt"
    threads:
        4
    shell:
        """
        samtools markdup -s -f {params.stats} -@ {threads} {input} {output} &> {log}
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
    # ValidateSamFile below, which only checks that the BAM is well-formed, not how well
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
        "bam.table"
    params:
        ref=config['reference']['key']
    run:
        with open(output[0],'w') as file:
            for sample in SAMPLES:
                file.write(f"{sample}\tbam_input/final/{sample}/{sample}.{params.ref}.bam\n")
