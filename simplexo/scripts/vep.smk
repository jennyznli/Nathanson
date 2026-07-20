from datetime import datetime
import os

# Get configuration values
PROJECT = config['project']['name']
RUN_DIR = config['project'].get('run_dir', '.')
CHROMOSOMES = list(range(1, 23))  # Chromosomes 1-22
date = datetime.now().strftime("%Y%m%d")

# Create run-specific directories
os.makedirs(f"{RUN_DIR}/logs", exist_ok=True)
os.makedirs(f"{RUN_DIR}/preprocess", exist_ok=True)
os.makedirs(f"{RUN_DIR}/output", exist_ok=True)
os.makedirs("reports", exist_ok=True)

wildcard_constraints:
    CHR='[1-9]|1[0-9]|2[0-2]'

rule run_vep_annotation:
    input:
        # Per-chromosome VEP reports
        expand(f"{RUN_DIR}/output/{{PROJECT}}.chr{{CHR}}.vep.report.csv", 
               PROJECT=PROJECT, CHR=CHROMOSOMES),
        # Final QC report
        f"reports/{PROJECT}_variant_qc_report.html"

# ========== VEP ANNOTATION PIPELINE ==========

rule no_sample:
    input:
        vcf=lambda wildcards: config['input']['qc_vcf'].format(CHR=wildcards.CHR)
    output:
        f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.no_sample.vcf"
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.no_sample.log"
    shell:
        """
        (bcftools view -G -Ov -o {output} {input.vcf})
        """

rule annotate_variants:
    input:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.no_sample.vcf"
    output:
        f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.no_sample.vep.vcf"
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.annotate_variants.log"
    shell:
        """
        (singularity run -H $PWD:/home \
        --bind /home/jennyzli/resources:/opt/vep/resources \
        --bind /home/jennyzli/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i {input.vcf} \
        -o {output} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta /opt/vep/resources/Genomes/Human/hg38/fa/Homo_sapiens_assembly38.fasta \
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
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.no_sample.vep.vcf"
    output:
        f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.vep.report.csv"
    params:
        blacklist=config['input']['blacklist']
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.parse_vep.log"
    shell:
        """
        (python vep_vcf_parser.py \
        -i {input.vcf} \
        -o {output} \
        -b {params.blacklist} \
        -m no_sample) 2>&1 | tee {log}
        """

rule qc_variants:
    input:
        reports=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.vep.report.csv", 
                      CHR=CHROMOSOMES)
    output:
        html=f"reports/{PROJECT}_variant_qc_report.html"
    params:
        input_pattern=f"{RUN_DIR}/output/{PROJECT}.chr*.vep.report.csv"
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.qc_variants.log"
    shell:
        """
        (python qc_variants.py {params.input_pattern} --json) 2>&1 | tee {log}
        """
