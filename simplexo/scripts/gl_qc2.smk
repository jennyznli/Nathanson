from datetime import datetime
import os

PROJECT = config['project']['name']
RUN_DIR = config['project'].get('run_dir', '.')
CHROMOSOMES = list(range(1, 23))  # Chromosomes 1-22
date = datetime.now().strftime("%Y%m%d")

# Create all necessary directories
os.makedirs(f"{RUN_DIR}/logs/chr_qc", exist_ok=True)
os.makedirs(f"{RUN_DIR}/vcf_stats", exist_ok=True)
os.makedirs(f"{RUN_DIR}/stats", exist_ok=True)  # Added this

wildcard_constraints:
    CHR='[0-9]+'  # More flexible constraint

rule all:
    input:
        expand(f"{RUN_DIR}/stats/chr{{CHR}}_gl.tags.qc.vep.report.csv", CHR=CHROMOSOMES),  # Fixed path
        f"{RUN_DIR}/vcf_stats/combined_qc_stats.txt",  # Fixed path
        f"{RUN_DIR}/stats/combined_vep_report.csv",  # Fixed path
        f"{RUN_DIR}/vcf_stats/qc_report.html",  # Fixed path
        f"{RUN_DIR}/variant_qc_report.html"  # Fixed path

rule vcf_stats:
    input:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.bcf"
    output:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.stats"
    threads: 1
    log:
        f"{RUN_DIR}/logs/chr_qc/chr{{CHR}}_stats.log"
    shell:
        """
        (bcftools stats -f PASS -s - {input} > {output}) 2>&1 | tee {log}
        """

rule no_sample:
    input:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.bcf"
    output:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.no_sample.vcf"
    threads: 1
    log:
        f"{RUN_DIR}/logs/chr_qc/chr{{CHR}}_no_sample.log"
    shell:
        """
        (bcftools view -G -Ov -o {output} {input}) 2>&1 | tee {log}
        """

rule annotate_variants:
    input:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.no_sample.vcf"
    output:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.no_sample.vep.vcf"
    params:
        vep_container=config['vep']['container'],
        vep_resources=config['vep']['resources_dir'],
        vep_cache=config['vep']['cache_dir'],
        fasta=config['vep']['fasta']
    threads: 1
    log:
        f"{RUN_DIR}/logs/chr_qc/chr{{CHR}}_vep.log"
    shell:
        """
        (singularity run -H $PWD:/home \
        --bind {params.vep_resources}:/opt/vep/resources \
        --bind {params.vep_cache}:/opt/vep/.vep \
        {params.vep_container} vep \
        --dir /opt/vep/.vep \
        -i {input} \
        -o {output} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta {params.fasta} \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz) 2>&1 | tee {log}
        """

rule parse_vep:
    input:
        f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.no_sample.vep.vcf"
    output:
        f"{RUN_DIR}/stats/chr{{CHR}}_gl.tags.qc.vep.report.csv"
    params:
        blacklist=config['input'].get('blacklist', 'blacklist.txt')
    threads: 1
    log:
        f"{RUN_DIR}/logs/chr_qc/chr{{CHR}}_parse.log"
    shell:
        """
        (python vep_vcf_parser.py -i {input} -o {output} -b {params.blacklist} -m no_sample) 2>&1 | tee {log}
        """

rule combine_stats:
    input:
        expand(f"{RUN_DIR}/chr{{CHR}}_gl.tags.qc.stats", CHR=CHROMOSOMES)
    output:
        f"{RUN_DIR}/vcf_stats/combined_qc_stats.txt"
    threads: 1
    log:
        f"{RUN_DIR}/logs/combine_stats.log"
    shell:
        """
        (cat {input} > {output}) 2>&1 | tee {log}
        """

rule qc_stats:
    input:
        f"{RUN_DIR}/vcf_stats/combined_qc_stats.txt"
    output:
         f"{RUN_DIR}/vcf_stats/qc_stats.csv",
         f"{RUN_DIR}/vcf_stats/qc_report.html"
    params:
        # metadata=config['input'].get('metadata', 'meta4qc.csv'),
        # bed=config['input'].get('bed', 'xgen_plus_spikein.Covered.slop10.hg38.bed')
    threads: 1
    log:
        f"{RUN_DIR}/logs/qc_stats.log"
    shell:
        """
        (python qc_stats.py --metadata {params.metadata} --bed {params.bed} {input}) 2>&1 | tee {log}
        """

rule combine_vep_reports:
    input:
        expand(f"{RUN_DIR}/stats/chr{{CHR}}_gl.tags.qc.vep.report.csv", CHR=CHROMOSOMES)
    output:
        f"{RUN_DIR}/stats/combined_vep_report.csv"
    threads: 1
    log:
        f"{RUN_DIR}/logs/combine_vep.log"
    run:
        import pandas as pd
        dfs = []
        for f in input:
            dfs.append(pd.read_csv(f))
        combined = pd.concat(dfs, ignore_index=True)
        combined.to_csv(output[0], index=False)

rule qc_variants:
    input:
        f"{RUN_DIR}/stats/combined_vep_report.csv"
    output:
        f"{RUN_DIR}/variant_qc_report.html"
    log:
        f"{RUN_DIR}/logs/qc_variants.log"
    shell:
        """
        (python qc_variants.py {input} --json) 2>&1 | tee {log}
        """