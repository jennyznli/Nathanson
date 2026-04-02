from datetime import datetime
import os

PROJECT = config['project']['name']
RUN_DIR = config['project'].get('run_dir', '.')
date = datetime.now().strftime("%Y%m%d")

GENE_INFO = {}
for gene, region_str in config['input']['gene_regions'].items():
    chrom, coords = region_str.split(':')
    chrom = chrom.replace('chr', '') 
    start, end = coords.replace(',', '').split('-') 
    GENE_INFO[gene] = {
        'chr': chrom,
        'start': int(start),
        'end': int(end),
        'region': region_str
    }

GENES = list(GENE_INFO.keys())
CHROMOSOMES = sorted(list(set([info['chr'] for info in GENE_INFO.values()])))

os.makedirs(f"{RUN_DIR}/logs", exist_ok=True)
os.makedirs(f"{RUN_DIR}/preprocess", exist_ok=True)
os.makedirs(f"{RUN_DIR}/output", exist_ok=True)
os.makedirs("reports", exist_ok=True)

wildcard_constraints:
    CHR='[XY0-9]+',
    GENE='[A-Za-z0-9_-]+'

def get_gene_chr(wildcards):
    return GENE_INFO[wildcards.GENE]['chr']

def get_gene_region(wildcards):
    return GENE_INFO[wildcards.GENE]['region']

rule all:
    input:
        [f"{RUN_DIR}/output/{PROJECT}.{gene}.vep.report.csv"
         for gene in GENES],
        expand(
            f"reports/{PROJECT}_{{GENE}}_variant_qc_report.html",
            PROJECT=PROJECT,
            GENE=GENES
        )

rule extract_gene_region:
    """Extract gene region (compressed) for VEP"""
    input:
        vcf=lambda wildcards: config['input']['qc_vcf'].format(
            CHR=get_gene_chr(wildcards)
        )
    output:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.extracted.vcf.gz"
    params:
        region=get_gene_region
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.extract_region.log"
    shell:
        """
        bcftools view -r {params.region} {input.vcf} -Oz -o {output.vcf} \
        2>&1 | tee {log}
        """

rule annotate_variants:
    """Run VEP on the extracted gene region"""
    input:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.extracted.vcf.gz"
    output:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.vep.vcf"
    threads: 2
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.annotate_variants.log"
    shell:
        """
        (singularity run -H $PWD:/home \
        --bind /home/jennyzli/resources:/opt/vep/resources \
        --bind /home/jennyzli/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i {input.vcf} \
        -o {output.vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --mane_select --canonical \
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
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz \
        ) 2>&1 | tee {log}

        rm -f {input.vcf}
        """

rule parse_vep:
    """Parse VEP output to CSV format"""
    input:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.vep.vcf"
    output:
        csv=f"{RUN_DIR}/output/{PROJECT}.{{GENE}}.vep.report.csv"
    params:
        blacklist=config['input']['blacklist']
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.parse_vep.log"
    shell:
        """
        export LD_LIBRARY_PATH=$HOME/software/htslib-1.21/lib:$LD_LIBRARY_PATH
        (python vep_vcf_parser.py \
        -i {input.vcf} \
        -o {output.csv} \
        -b {params.blacklist} \
        -m cohort) 2>&1 | tee {log}
        """

rule qc_variants:
    """Generate QC report for each gene"""
    input:
        report=f"{RUN_DIR}/output/{PROJECT}.{{GENE}}.vep.report.csv"
    output:
        html=f"reports/{PROJECT}_{{GENE}}_variant_qc_report.html"
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.qc_variants.log"
    shell:
        """
        (python qc_variants.py {input.report} --json) 2>&1 | tee {log}
        """
