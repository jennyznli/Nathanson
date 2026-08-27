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
os.makedirs(f"{RUN_DIR}/figures", exist_ok=True)
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
            GENE=GENES
        ),
        expand(
            f"{RUN_DIR}/output/{PROJECT}.{{GENE}}.vep.report_{{GENE}}_qc.csv",
            GENE=GENES
        ),
        expand(
            f"{RUN_DIR}/figures/{PROJECT}.{{GENE}}.vep.report_{{GENE}}_qc_histograms.pdf",
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
        bcftools view -f PASS -r {params.region} {input.vcf} -Oz -o {output.vcf} \
        2>&1 | tee {log}
        """

rule annotate_variants:
    """Run VEP on the extracted gene region"""
    input:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.extracted.vcf.gz"
    output:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.vep.vcf"
    threads: 16
    resources:
        mem_mb=48000
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
        --fork {threads} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --mane --canonical \
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
        """

rule filter_annotated_pass:
    """Filter the ALREADY-annotated VEP VCF down to PASS-only sites.
    Operates on the existing preprocess/*.vep.vcf so annotate_variants
    (slow, container/VEP) never has to re-run for this pass -- only used
    to salvage a run that was extracted/annotated before -f PASS was
    added to extract_gene_region."""
    input:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.vep.vcf"
    output:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.vep.pass.vcf"
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.filter_annotated_pass.log"
    shell:
        """
        bcftools view -f PASS {input.vcf} -Ov -o {output.vcf} \
        2>&1 | tee {log}
        """

rule parse_vep:
    """Parse VEP output to CSV format"""
    input:
        vcf=f"{RUN_DIR}/preprocess/{PROJECT}.{{GENE}}.vep.pass.vcf"
    output:
        csv=f"{RUN_DIR}/output/{PROJECT}.{{GENE}}.vep.report.csv"
    resources:
        mem_mb=32000
    params:
        blacklist=config['input']['blacklist']
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.parse_vep.log"
    shell:
        """
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vep_parser
        python scripts/vep_vcf_parser.py \
        -i {input.vcf} \
        -o {output.csv} \
        -b {params.blacklist} \
        -m cohort 2>&1 | tee {log}
        """

rule filter_germline_qc:
    """QC-filter a gene's VEP report to high-confidence germline calls
    (R/germline_qc.R::qc_germline_gene) and write a diagnostic histogram PDF"""
    input:
        csv=f"{RUN_DIR}/output/{PROJECT}.{{GENE}}.vep.report.csv"
    output:
        csv=f"{RUN_DIR}/output/{PROJECT}.{{GENE}}.vep.report_{{GENE}}_qc.csv",
        pdf=f"{RUN_DIR}/figures/{PROJECT}.{{GENE}}.vep.report_{{GENE}}_qc_histograms.pdf"
    params:
        gene=lambda wildcards: wildcards.GENE,
        log_dir=f"{RUN_DIR}/logs"
    resources:
        mem_mb=4000
    log:
        f"{RUN_DIR}/logs/{PROJECT}.{{GENE}}.filter_germline_qc.log"
    shell:
        """
        Rscript scripts/run_germline_qc_cli.R \
        --input {input.csv} \
        --gene {params.gene} \
        --output_csv {output.csv} \
        --output_pdf {output.pdf} \
        --log_dir {params.log_dir} \
        2>&1 | tee {log}
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
        (python scripts/qc_variants.py {input.report} --json) 2>&1 | tee {log}
        """
