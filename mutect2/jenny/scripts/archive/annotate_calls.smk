"""
Post-calling annotation for the tumor Mutect2 cohort (deep targeted panel, no chromosome
scatter -- the whole point of a targeted assay is that per-sample call volume is small
enough to normalize/annotate/parse in one shot per sample).

Input is each sample's already-filtered Mutect2 VCF, not a genotyped BCF from an
upstream array/WES pipeline.

"""

import os

os.makedirs("results/logs/annotate", exist_ok=True)
os.makedirs("results/annotate/norm", exist_ok=True)
os.makedirs("results/annotate/vep", exist_ok=True)
os.makedirs("results/annotate/csv", exist_ok=True)

REF_FASTA = config["reference"]["fasta"]
TARGETS_BED = config["resources"]["targets_bed"]


def _load_sample_tsv(path):
    samples = {}
    if not os.path.exists(path):
        return samples
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            samples[parts[0].strip()] = parts[1].strip()
    return samples


# Scope: the tumor calling cohort only (not the PoN normals)
TUMOR_SAMPLES = _load_sample_tsv(config["samples"]["tumor_tsv"])
if not TUMOR_SAMPLES:
    raise ValueError(f"No samples found in {config['samples']['tumor_tsv']}")

wildcard_constraints:
    sample="|".join(__import__("re").escape(s) for s in TUMOR_SAMPLES)


rule annotate_all:
    input:
        expand("results/annotate/csv/{sample}.vep.csv", sample=TUMOR_SAMPLES),
        "results/annotate/all_samples.vep.csv"


rule normalize_and_set_id:
    input:
        vcf="results/calls/{sample}/{sample}.filtered.vcf.gz"
    output:
        vcf="results/annotate/norm/{sample}.norm.vcf.gz",
        csi="results/annotate/norm/{sample}.norm.vcf.gz.csi"
    log:
        "results/logs/annotate/{sample}.normalize.log"
    shell:
        r"""
        exec > {log} 2>&1
        bcftools norm -m-any -f {REF_FASTA} {input.vcf} \
            | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' \
            | bcftools view -Oz -o {output.vcf}
        bcftools index -f {output.vcf}
        """


rule vep_annotate:
    input:
        "results/annotate/norm/{sample}.norm.vcf.gz"
    output:
        "results/annotate/vep/{sample}.vep.vcf"
    threads: 8
    resources:
        mem_mb=24000
    log:
        "results/logs/annotate/{sample}.vep.log"
    shell:
        r"""
        (export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache

        mkdir -p /scratch/$USER/bin/{wildcards.sample}
        cp /appl/samtools-1.21/bin/samtools /scratch/$USER/bin/{wildcards.sample}/samtools
        chmod +x /scratch/$USER/bin/{wildcards.sample}/samtools

        LOFTEE_REPO="/home/jennyzli/.vep/Plugins/loftee"
        LOFTEE_DIR="${{LOFTEE_REPO}}/GRCh38"

        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        --bind /home/jennyzli/resources:/opt/vep/resources \
        --bind /home/jennyzli/.vep:/opt/vep/.vep \
        --bind /scratch/$USER/bin/{wildcards.sample}/samtools:/usr/local/bin/samtools \
        --bind ${{LOFTEE_REPO}}:${{LOFTEE_REPO}} \
        --bind ${{LOFTEE_REPO}}/LoF.pm:/plugins/LoF.pm \
        --bind ${{LOFTEE_REPO}}/ancestral.pm:/plugins/ancestral.pm \
        --bind ${{LOFTEE_REPO}}/context.pm:/plugins/context.pm \
        --bind ${{LOFTEE_REPO}}/de_novo_donor.pl:/plugins/de_novo_donor.pl \
        --bind ${{LOFTEE_REPO}}/extended_splice.pl:/plugins/extended_splice.pl \
        --bind ${{LOFTEE_REPO}}/gerp_dist.pl:/plugins/gerp_dist.pl \
        --bind ${{LOFTEE_REPO}}/loftee_splice_utils.pl:/plugins/loftee_splice_utils.pl \
        --bind ${{LOFTEE_REPO}}/splice_site_scan.pl:/plugins/splice_site_scan.pl \
        --bind ${{LOFTEE_REPO}}/svm.pl:/plugins/svm.pl \
        --bind ${{LOFTEE_REPO}}/TissueExpression.pm:/plugins/TissueExpression.pm \
        --bind ${{LOFTEE_REPO}}/utr_splice.pl:/plugins/utr_splice.pl \
        --bind ${{LOFTEE_REPO}}/maxEntScan:/plugins/maxEntScan \
        /appl/containers/ensembl-vep_release_116.0.sif vep \
        --dir /opt/vep/.vep \
        -i $PWD/{input} \
        -o $PWD/{output} \
        --fork {threads} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical --mane \
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
        --plugin LoF,loftee_path:${{LOFTEE_REPO}},human_ancestor_fa:${{LOFTEE_DIR}}/human_ancestor.fa.gz,conservation_file:${{LOFTEE_DIR}}/loftee.sql,gerp_bigwig:${{LOFTEE_DIR}}/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
        --cache_version 116) 2>&1 | tee {log}
        """


rule parse_vep_csv:
    input:
        "results/annotate/vep/{sample}.vep.vcf"
    output:
        "results/annotate/csv/{sample}.vep.csv"
    log:
        "results/logs/annotate/{sample}.parse.log"
    shell:
        r"""
        (source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vep_parser
        python scripts/vep_vcf_parser.py -i {input} -o {output} -m cohort) 2>&1 | tee {log}
        """


rule merge_samples:
    input:
        expand("results/annotate/csv/{sample}.vep.csv", sample=TUMOR_SAMPLES)
    output:
        "results/annotate/all_samples.vep.csv"
    log:
        "results/logs/annotate/merge_samples.log"
    shell:
        r"""
        exec > {log} 2>&1
        awk 'FNR==1 && NR!=1{{next}}{{print}}' {input} > {output}
        """
