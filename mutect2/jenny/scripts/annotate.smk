"""
Caller-agnostic post-calling annotation for the tumor cohort. 

    snakemake -s scripts/annotate.smk --config caller=mutect2
    snakemake -s scripts/annotate.smk --config caller=vardict
    snakemake -s scripts/annotate.smk --config caller=vardict outdir=/some/other/path

"""

import os

DEFAULT_OUTDIR = {
    "mutect2": "mutect2",
    "vardict": "vardict",
}

CALLER = config.get("caller")
if CALLER not in DEFAULT_OUTDIR:
    raise ValueError(f"config['caller'] must be one of {list(DEFAULT_OUTDIR)}, got {CALLER!r}")

OUTDIR = config.get("outdir", DEFAULT_OUTDIR[CALLER])

# if it's mutect2, get the right subanalysis (pon vs no pon)
if CALLER == "mutect2":
    ANALYSIS_LABEL = (config.get("analysis") or {}).get("label")
    if not ANALYSIS_LABEL:
        raise ValueError(
            "config['analysis']['label'] is required for caller=mutect2 -- "
            "mutect2.smk writes filtered VCFs to mutect2/calls/{label}/"
        )
    BASENAME = f"{{sample}}.{ANALYSIS_LABEL}"
    RAW_VCF = f"{OUTDIR}/calls/{ANALYSIS_LABEL}/{{sample}}/{BASENAME}.filtered.vcf.gz"
    LOGS_DIR = f"{OUTDIR}/logs/{ANALYSIS_LABEL}"
    NORM_DIR = f"{OUTDIR}/norm/{ANALYSIS_LABEL}"
    VEP_DIR = f"{OUTDIR}/vep/{ANALYSIS_LABEL}"
    CSV_DIR = f"{OUTDIR}/csv/{ANALYSIS_LABEL}"
    MERGED_CSV = f"{OUTDIR}/all_samples.{CALLER}.{ANALYSIS_LABEL}.vep.csv"
    ONCOKB_CSV = f"{OUTDIR}/all_samples.{CALLER}.{ANALYSIS_LABEL}.vep.oncokb.csv"
    SCRATCH_TAG = f"{CALLER}/{ANALYSIS_LABEL}"
else:
    BASENAME = "{sample}"
    RAW_VCF = f"{OUTDIR}/calls/{{sample}}/{BASENAME}.vcf.gz"
    LOGS_DIR = f"{OUTDIR}/logs"
    NORM_DIR = f"{OUTDIR}/norm"
    VEP_DIR = f"{OUTDIR}/vep"
    CSV_DIR = f"{OUTDIR}/csv"
    MERGED_CSV = f"{OUTDIR}/all_samples.{CALLER}.vep.csv"
    ONCOKB_CSV = f"{OUTDIR}/all_samples.{CALLER}.vep.oncokb.csv"
    SCRATCH_TAG = CALLER

os.makedirs(LOGS_DIR, exist_ok=True)
os.makedirs(NORM_DIR, exist_ok=True)
os.makedirs(VEP_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)

# OUTDIR is absolute for callers like vardict (default outside $PWD) -- in that
# case don't prepend $PWD to already-absolute vep_annotate paths, and bind
# OUTDIR itself into the container since only $PWD is bound by default.
PWD_PREFIX = "" if os.path.isabs(OUTDIR) else "$PWD/"
OUTDIR_BIND = f"--bind {OUTDIR}:{OUTDIR} \\\n        " if os.path.isabs(OUTDIR) else ""

REF_FASTA = config["reference"]["fasta"]
TARGETS_BED = config["resources"]["targets_bed"]

# one OncoKB cache shared by every caller/label run from this dir; token comes from
# $ONCOKB_TOKEN or ~/.oncokb_token, never config
ONCOKB_CACHE = (config.get("oncokb") or {}).get("cache", "oncokb/cache.GRCh38.tsv")


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


# tumor calling cohort only (not the PoN normals)
TUMOR_SAMPLES = _load_sample_tsv(config["samples"]["cases_tsv"])
if not TUMOR_SAMPLES:
    raise ValueError(f"No samples found in {config['samples']['cases_tsv']}")

wildcard_constraints:
    sample="|".join(__import__("re").escape(s) for s in TUMOR_SAMPLES)


# oncokb_annotate needs internet -- run it where snakemake runs instead of submitting to LSF
localrules: annotate_all, oncokb_annotate


rule annotate_all:
    input:
        expand(f"{CSV_DIR}/{BASENAME}.vep.csv", sample=TUMOR_SAMPLES),
        MERGED_CSV,
        ONCOKB_CSV


rule normalize_and_set_id:
    input:
        vcf=RAW_VCF
    output:
        vcf=f"{NORM_DIR}/{BASENAME}.norm.vcf.gz",
        csi=f"{NORM_DIR}/{BASENAME}.norm.vcf.gz.csi"
    log:
        f"{LOGS_DIR}/{{sample}}.normalize.log"
    shell:
        r"""
        exec > {log} 2>&1
        bcftools norm -m -any -f {REF_FASTA} {input.vcf} \
            | bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' \
            | bcftools view -Oz -o {output.vcf}
        bcftools index -f {output.vcf}
        """


rule vep_annotate:
    input:
        f"{NORM_DIR}/{BASENAME}.norm.vcf.gz"
    output:
        f"{VEP_DIR}/{BASENAME}.vep.vcf"
    threads: 8
    resources:
        mem_mb=24000
    log:
        f"{LOGS_DIR}/{{sample}}.vep.log"
    shell:
        r"""
        (export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache

        mkdir -p /scratch/$USER/bin/{SCRATCH_TAG}/{wildcards.sample}
        cp /appl/samtools-1.21/bin/samtools /scratch/$USER/bin/{SCRATCH_TAG}/{wildcards.sample}/samtools
        chmod +x /scratch/$USER/bin/{SCRATCH_TAG}/{wildcards.sample}/samtools

        LOFTEE_REPO="/home/jennyzli/.vep/Plugins/loftee"
        LOFTEE_DIR="${{LOFTEE_REPO}}/GRCh38"

        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        {OUTDIR_BIND}--bind /home/jennyzli/resources:/opt/vep/resources \
        --bind /home/jennyzli/.vep:/opt/vep/.vep \
        --bind /scratch/$USER/bin/{SCRATCH_TAG}/{wildcards.sample}/samtools:/usr/local/bin/samtools \
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
        -i {PWD_PREFIX}{input} \
        -o {PWD_PREFIX}{output} \
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
        f"{VEP_DIR}/{BASENAME}.vep.vcf"
    output:
        f"{CSV_DIR}/{BASENAME}.vep.csv"
    log:
        f"{LOGS_DIR}/{{sample}}.parse.log"
    shell:
        r"""
        (source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vep_parser
        python scripts/vep_vcf_parser.py -i {input} -o {output} -m cohort) 2>&1 | tee {log}
        """


rule merge_samples:
    input:
        expand(f"{CSV_DIR}/{BASENAME}.vep.csv", sample=TUMOR_SAMPLES)
    output:
        MERGED_CSV
    log:
        f"{LOGS_DIR}/merge_samples.log"
    shell:
        r"""
        exec > {log} 2>&1
        awk 'FNR==1 && NR!=1{{next}}{{print}}' {input} > {output}
        """


rule oncokb_annotate:
    input:
        MERGED_CSV
    output:
        ONCOKB_CSV
    log:
        f"{LOGS_DIR}/oncokb_annotate.log"
    shell:
        r"""
        (source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vep_parser
        python scripts/oncokb_annotate.py -i {input} -o {output} --cache {ONCOKB_CACHE}) 2>&1 | tee {log}
        """
