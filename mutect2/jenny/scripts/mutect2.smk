"""
Tumor-only Mutect2 calling, no scatter
output: one VCF and index w filtering applied
"""

import os
import re

# ---------------------------------------------------------------------------
# Config / resources
# ---------------------------------------------------------------------------

REF_FASTA = config["reference"]["fasta"]
REF_FAI = config["reference"]["fai"]
REF_DICT = config["reference"]["dict"]

GNOMAD = config["resources"]["gnomad"]
VARIANTS_FOR_CONTAMINATION = config["resources"]["variants_for_contamination"]
TARGETS_BED = config["resources"]["targets_bed"]
MASKED_INTERVALS = config["resources"].get("masked_intervals", "") or ""

PON_CFG = config["pon"]
PARAMS = config["parameters"]

USE_PON = bool(PON_CFG.get("use_pon", True))
EXISTING_PON = (PON_CFG.get("existing_pon_vcf") or "").strip()
BUILD_PON = USE_PON and not bool(EXISTING_PON)
PON_VCF = EXISTING_PON if EXISTING_PON else f"mutect2/pon/{PON_CFG['name']}.vcf.gz"

MASKED_ARG = f"-XL {MASKED_INTERVALS}" if MASKED_INTERVALS else ""

ANALYSIS_CFG = config.get("analysis") or {}
ANALYSIS_LABEL = ANALYSIS_CFG.get("label")
if not ANALYSIS_LABEL:
    raise ValueError(
        "config['analysis']['label'] is required (e.g. 'with_pon' / 'without_pon') -- "
        "keeps different with/without-PoN comparison runs from overwriting each other's "
        "mutect2/calls/ output"
    )

TUMOR_DIR = f"mutect2/calls/{ANALYSIS_LABEL}"
PON_NORMAL_DIR = "mutect2/pon/normals"

# ---------------------------------------------------------------------------
# Sample sheets
# ---------------------------------------------------------------------------

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
            sample_id, path_ = parts[0].strip(), parts[1].strip()
            samples[sample_id] = path_
    return samples

NORMAL_SAMPLES = _load_sample_tsv(config["samples"]["normals_tsv"]) if BUILD_PON else {}
TUMOR_SAMPLES = _load_sample_tsv(config["samples"]["cases_tsv"])

_overlap = set(NORMAL_SAMPLES) & set(TUMOR_SAMPLES)
if _overlap:
    raise ValueError(f"Sample IDs appear in both normals_tsv and cases_tsv: {sorted(_overlap)}")

SAMPLE_PATH = {**NORMAL_SAMPLES, **TUMOR_SAMPLES}

if BUILD_PON and not NORMAL_SAMPLES:
    raise ValueError(
        "PoN construction requested, but no normal samples were loaded "
        "from samples.normals_tsv"
    )

if not TUMOR_SAMPLES:
    raise ValueError(
        "No tumor samples were loaded from samples.cases_tsv"
    )

wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLE_PATH),


NORMAL_SAMPLE_PATTERN = "|".join(re.escape(s) for s in NORMAL_SAMPLES) if NORMAL_SAMPLES else "NOMATCH"
TUMOR_SAMPLE_PATTERN = "|".join(re.escape(s) for s in TUMOR_SAMPLES) if TUMOR_SAMPLES else "NOMATCH"


def sample_path(wildcards):
    return SAMPLE_PATH[wildcards.sample]


def index_for(wildcards):
    p = SAMPLE_PATH[wildcards.sample]
    if p.endswith(".cram"):
        for cand in (p + ".crai", re.sub(r"\.cram$", ".crai", p)):
            if os.path.exists(cand):
                return cand
        return p + ".crai"
    if p.endswith(".bam"):
        for cand in (p + ".bai", re.sub(r"\.bam$", ".bai", p)):
            if os.path.exists(cand):
                return cand
        return p + ".bai"
    raise ValueError(f"Unrecognized read file extension for sample path: {p}")


# ---------------------------------------------------------------------------
# Top-level target
# ---------------------------------------------------------------------------

def all_targets(wildcards):
    targets = [f"{TUMOR_DIR}/{s}/{s}.{ANALYSIS_LABEL}.filtered.vcf.gz" for s in TUMOR_SAMPLES]
    if BUILD_PON:
        targets.append(PON_VCF)
    return targets


rule all:
    input:
        all_targets


# ---------------------------------------------------------------------------
# PoN-normal calling + filtering 
# ---------------------------------------------------------------------------

rule mutect2_pon_normal:
    wildcard_constraints:
        sample=NORMAL_SAMPLE_PATTERN,
    input:
        ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
        reads=sample_path,
        index=index_for,
        targets=TARGETS_BED,
    output:
        vcf=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.unfiltered.vcf.gz",
        idx=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.unfiltered.vcf.gz.tbi",
        stats=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.unfiltered.vcf.gz.stats",
    log:
        "mutect2/logs/mutect2_pon_normal/{sample}.log",
    params:
        extra=PARAMS.get("m2_extra_args", "") or "",
        masked=MASKED_ARG,
    threads: 4
    resources:
        mem_mb=6000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx5g" Mutect2 \
            -R {input.ref} \
            -I {input.reads} \
            -L {input.targets} \
            {params.masked} \
            --max-mnp-distance 0 \
            {params.extra} \
            -O {output.vcf}
        """


rule filter_mutect_calls_pon_normal:
    wildcard_constraints:
        sample=NORMAL_SAMPLE_PATTERN,
    input:
        vcf=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.unfiltered.vcf.gz",
        ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
        stats=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.unfiltered.vcf.gz.stats",
    output:
        vcf=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.filtered.vcf.gz",
        idx=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.filtered.vcf.gz.tbi",
        stats=f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.filtering.stats",
    log:
        "mutect2/logs/filter_mutect_calls_pon_normal/{sample}.log",
    params:
        extra=PARAMS.get("m2_extra_filtering_args", "") or "",
    resources:
        mem_mb=5000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx4g" FilterMutectCalls \
            -V {input.vcf} \
            -R {input.ref} \
            -O {output.vcf} \
            -stats {input.stats} \
            --filtering-stats {output.stats} {params.extra}
        """


# ---------------------------------------------------------------------------
# Tumor-cohort calling 
# ---------------------------------------------------------------------------

def pon_arg(wildcards):
    return f"-pon {PON_VCF}" if USE_PON else ""


rule mutect2_tumor:
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
        reads=sample_path,
        index=index_for,
        targets=TARGETS_BED,
        pon=[PON_VCF, PON_VCF + ".tbi"] if USE_PON else [],
    output:
        vcf=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.unfiltered.vcf.gz",
        idx=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.unfiltered.vcf.gz.tbi",
        stats=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.unfiltered.vcf.gz.stats",
        f1r2=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.f1r2.tar.gz",
    log:
        f"mutect2/logs/{ANALYSIS_LABEL}/mutect2_tumor/{{sample}}.log",
    params:
        pon_arg=pon_arg,
        extra=PARAMS.get("m2_extra_args", "") or "",
        masked=MASKED_ARG,
    threads: 4
    resources:
        mem_mb=6000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx5g" Mutect2 \
            -R {input.ref} \
            -I {input.reads} \
            -L {input.targets} \
            {params.masked} \
            {params.pon_arg} \
            --germline-resource {GNOMAD} \
            --f1r2-tar-gz {output.f1r2} \
            {params.extra} \
            -O {output.vcf}
        """


# tabulates AF pileup data at known common SNP sites, used later on by CalculateContamination
rule get_pileup_summaries:
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        ref=REF_FASTA,
        reads=sample_path,
        index=index_for,
        targets=TARGETS_BED,
        variants_for_contamination=VARIANTS_FOR_CONTAMINATION,
    output:
        table=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.pileups.table",
    log:
        f"mutect2/logs/{ANALYSIS_LABEL}/get_pileup_summaries/{{sample}}.log",
    resources:
        mem_mb=4000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx3g" GetPileupSummaries \
            -R {input.ref} \
            -I {input.reads} \
            --interval-set-rule INTERSECTION \
            -L {input.targets} \
            -V {input.variants_for_contamination} \
            -L {input.variants_for_contamination} \
            -O {output.table}
        """


# learning step of orientation bias mixture model
rule learn_read_orientation_model:
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        f1r2=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.f1r2.tar.gz",
    output:
        f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.artifact-priors.tar.gz",
    log:
        f"mutect2/logs/{ANALYSIS_LABEL}/learn_read_orientation_model/{{sample}}.log",
    resources:
        mem_mb=7000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx6g" LearnReadOrientationModel -I {input.f1r2} -O {output}
        """

# Tumor-only: no -matched normal pileups.
rule calculate_contamination:
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        tumor_pileups=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.pileups.table",
    output:
        contamination=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.contamination.table",
        segments=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.segments.table",
    log:
        f"mutect2/logs/{ANALYSIS_LABEL}/calculate_contamination/{{sample}}.log",
    resources:
        mem_mb=4000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx3g" CalculateContamination \
            -I {input.tumor_pileups} \
            -O {output.contamination} \
            --tumor-segmentation {output.segments}
        """


rule filter_mutect_calls_tumor:
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        vcf=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.unfiltered.vcf.gz",
        ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
        stats=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.unfiltered.vcf.gz.stats",
        priors=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.artifact-priors.tar.gz",
        contamination=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.contamination.table",
        segments=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.segments.table",
    output:
        vcf=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.filtered.vcf.gz",
        idx=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.filtered.vcf.gz.tbi",
        stats=f"{TUMOR_DIR}/{{sample}}/{{sample}}.{ANALYSIS_LABEL}.filtering.stats",
    log:
        f"mutect2/logs/{ANALYSIS_LABEL}/filter_mutect_calls_tumor/{{sample}}.log",
    params:
        extra=PARAMS.get("m2_extra_filtering_args", "") or "",
    resources:
        mem_mb=5000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx4g" FilterMutectCalls \
            -V {input.vcf} \
            -R {input.ref} \
            -O {output.vcf} \
            -stats {input.stats} \
            --filtering-stats {output.stats} \
            --contamination-table {input.contamination} \
            --tumor-segmentation {input.segments} \
            --ob-priors {input.priors} {params.extra}
        """


# ---------------------------------------------------------------------------
# Panel of Normals construction
# ---------------------------------------------------------------------------

if BUILD_PON:

    rule genomicsdb_and_panel:
        input:
            vcfs=expand(f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.filtered.vcf.gz", sample=list(NORMAL_SAMPLES)),
            tbis=expand(f"{PON_NORMAL_DIR}/{{sample}}/{{sample}}.filtered.vcf.gz.tbi", sample=list(NORMAL_SAMPLES)),
            ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
            gnomad=GNOMAD,
            targets=TARGETS_BED,
        output:
            vcf=PON_VCF,
            idx=PON_VCF + ".tbi",
        log:
            "mutect2/logs/genomicsdb_and_panel/genomicsdb_and_panel.log",
        params:
            vcf_args=lambda wc, input: " -V ".join(input.vcfs),
            masked=MASKED_ARG,
        resources:
            mem_mb=9000,
        shell:
            r"""
            exec > {log} 2>&1
            if [ -w /scratch ]; then SCR=/scratch; else SCR=/tmp; fi
            mkdir -p $SCR/$(id -un)
            WORK=$(mktemp -d $SCR/$(id -un)/gdb_panel.XXXXXX)
            trap 'rm -rf "$WORK"' EXIT
            mkdir -p $WORK/tmp
            hostname; df -h $SCR /tmp

            gatk --java-options "-Xmx8g -Djava.io.tmpdir=$WORK/tmp" GenomicsDBImport \
                --genomicsdb-workspace-path $WORK/db \
                --tmp-dir $WORK/tmp \
                -R {input.ref} \
                -V {params.vcf_args} \
                -L {input.targets} \
                {params.masked}

            gatk --java-options "-Xmx8g -Djava.io.tmpdir=$WORK/tmp" CreateSomaticPanelOfNormals \
                -R {input.ref} \
                --germline-resource {input.gnomad} \
                -V gendb://$WORK/db \
                -O {output.vcf}
            """
