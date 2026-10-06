"""
Tumor-only calling pipeline (Mutect2) w/ scattering
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

SCATTER_COUNT = config["scatter_count"]

PON_CFG = config["pon"]
TUMOR_CFG = config["tumor_calling"]

EXISTING_PON = (PON_CFG.get("existing_pon_vcf") or "").strip()
BUILD_PON = not bool(EXISTING_PON)
PON_VCF = EXISTING_PON if EXISTING_PON else f"results/pon/{PON_CFG['name']}.vcf.gz"

MASKED_ARG = f"-XL {MASKED_INTERVALS}" if MASKED_INTERVALS else ""


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
TUMOR_SAMPLES = _load_sample_tsv(config["samples"]["tumor_tsv"])

_overlap = set(NORMAL_SAMPLES) & set(TUMOR_SAMPLES)
if _overlap:
    raise ValueError(f"Sample IDs appear in both normals_tsv and tumor_tsv: {sorted(_overlap)}")

SAMPLE_META = {}
for s, p in NORMAL_SAMPLES.items():
    SAMPLE_META[s] = {"path": p, "group": "pon"}
for s, p in TUMOR_SAMPLES.items():
    SAMPLE_META[s] = {"path": p, "group": "tumor"}

if not SAMPLE_META:
    raise ValueError("No samples found -- fill in samples/normals.tsv and/or samples/tumor_samples.tsv")

wildcard_constraints:
    sample = "|".join(re.escape(s) for s in SAMPLE_META),


def sample_group(sample):
    return SAMPLE_META[sample]["group"]


def is_tumor(wildcards):
    return sample_group(wildcards.sample) == "tumor"


TUMOR_SAMPLE_PATTERN = "|".join(re.escape(s) for s in TUMOR_SAMPLES) if TUMOR_SAMPLES else "NOMATCH"


def sample_path(wildcards):
    return SAMPLE_META[wildcards.sample]["path"]


def index_for(wildcards):
    p = SAMPLE_META[wildcards.sample]["path"]
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
    targets = [f"results/calls/{s}/{s}.filtered.vcf.gz" for s in TUMOR_SAMPLES]
    if BUILD_PON:
        targets.append(PON_VCF)
    return targets


rule all:
    input:
        all_targets


# ---------------------------------------------------------------------------
# Interval scatter 
# ---------------------------------------------------------------------------

checkpoint split_intervals:
    input:
        ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT, targets=TARGETS_BED,
    output:
        directory("results/intervals"),
    log:
        "results/logs/split_intervals/split_intervals.log",
    params:
        scatter=SCATTER_COUNT,
        masked=MASKED_ARG,
    threads: 1
    resources:
        mem_mb=6000,
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {output}
        gatk --java-options "-Xmx5g" SplitIntervals \
            -R {input.ref} \
            -L {input.targets} \
            {params.masked} \
            -scatter {params.scatter} \
            -O {output}
        """


def shard_ids(_wildcards=None):
    ck = checkpoints.split_intervals.get()
    (shards,) = glob_wildcards(os.path.join(ck.output[0], "{shard}.interval_list"))
    return sorted(shards)


# ---------------------------------------------------------------------------
# Per-sample tumor-only calling + filtering (shared by pon and tumor cohorts)
# ---------------------------------------------------------------------------

rule get_sample_name:
    input:
        ref=REF_FASTA,
        reads=sample_path,
        index=index_for,
    output:
        "results/calls/{sample}/tumor_name.txt",
    log:
        "results/logs/get_sample_name/{sample}.log",
    resources:
        mem_mb=4000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx3g" GetSampleName -R {input.ref} -I {input.reads} -O {output} -encode
        """

def pon_inputs(wildcards):
    if sample_group(wildcards.sample) == "tumor":
        return [PON_VCF, PON_VCF + ".tbi"]
    return []


def pon_arg(wildcards):
    return f"-pon {PON_VCF}" if sample_group(wildcards.sample) == "tumor" else ""


def germline_resource_arg(wildcards):
    return f"--germline-resource {GNOMAD}" if sample_group(wildcards.sample) == "tumor" else ""


def max_pop_af_arg(wildcards):
    if sample_group(wildcards.sample) == "tumor":
        return f"--max-population-af {TUMOR_CFG['max_population_af']}"
    return ""


def max_mnp_distance(wildcards):
    # PoN construction hard-codes 0 (GenomicsDBImport can't handle MNPs), matching broad/mutect2_pon.wdl
    if sample_group(wildcards.sample) == "pon":
        return 0
    return TUMOR_CFG["max_mnp_distance"]


rule mutect2_scatter:
    input:
        ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
        reads=sample_path,
        index=index_for,
        tumor_name="results/calls/{sample}/tumor_name.txt",
        intervals="results/intervals/{shard}.interval_list",
        pon=pon_inputs,
    output:
        vcf="results/calls/{sample}/shards/{shard}.unfiltered.vcf.gz",
        idx="results/calls/{sample}/shards/{shard}.unfiltered.vcf.gz.tbi",
        stats="results/calls/{sample}/shards/{shard}.unfiltered.vcf.gz.stats",
        f1r2="results/calls/{sample}/shards/{shard}.f1r2.tar.gz",
    log:
        "results/logs/mutect2_scatter/{sample}/{shard}.log",
    params:
        pon_arg=pon_arg,
        germline_arg=germline_resource_arg,
        max_pop_af_arg=max_pop_af_arg,
        max_mnp=max_mnp_distance,
        min_af=TUMOR_CFG["minimum_allele_fraction"],
        lod_emit=TUMOR_CFG["tumor_lod_to_emit"],
        lod_call=TUMOR_CFG["initial_tumor_lod"],
        extra=TUMOR_CFG.get("m2_extra_args", "") or "",
    threads: 2
    resources:
        mem_mb=6000,
    shell:
        r"""
        exec > {log} 2>&1
        tumor_name=$(cat {input.tumor_name})
        gatk --java-options "-Xmx5g" Mutect2 \
            -R {input.ref} \
            -I {input.reads} \
            -tumor "$tumor_name" \
            -L {input.intervals} \
            {params.pon_arg} \
            {params.germline_arg} \
            {params.max_pop_af_arg} \
            --minimum-allele-fraction {params.min_af} \
            --tumor-lod-to-emit {params.lod_emit} \
            --initial-tumor-lod {params.lod_call} \
            --max-mnp-distance {params.max_mnp} \
            --native-pair-hmm-threads {threads} \
            --f1r2-tar-gz {output.f1r2} \
            {params.extra} \
            -O {output.vcf}
        """


rule get_pileup_summaries:
    # tumor cohort only -- see module docstring / broad/mutect2_pon.wdl: normals never
    # get variants_for_contamination wired in, so this step doesn't apply to them.
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        ref=REF_FASTA,
        reads=sample_path,
        index=index_for,
        intervals="results/intervals/{shard}.interval_list",
        resource=VARIANTS_FOR_CONTAMINATION,
    output:
        table="results/calls/{sample}/shards/{shard}.pileups.table",
    log:
        "results/logs/get_pileup_summaries/{sample}/{shard}.log",
    resources:
        mem_mb=4000,
    shell:
        # A shard with no sites intersecting the contamination resource can make
        # GetPileupSummaries fail/produce nothing (see broad/bick's WDL comments on
        # this same behavior) -- write a header-only table in that case so gathering
        # downstream doesn't break.
        r"""
        exec > {log} 2>&1
        set +e
        gatk --java-options "-Xmx3g" GetPileupSummaries \
            -R {input.ref} \
            -I {input.reads} \
            --interval-set-rule INTERSECTION \
            -L {input.intervals} \
            -V {input.resource} \
            -L {input.resource} \
            -O {output.table}
        rc=$?
        if [ ! -s {output.table} ]; then
            printf "contig\tposition\tref_count\talt_count\tother_alt_count\tallele_frequency\n" > {output.table}
        fi
        exit 0
        """


def _shard_expand(pattern):
    def _inner(wildcards):
        return expand(pattern, sample=wildcards.sample, shard=shard_ids())
    return _inner


rule merge_vcfs:
    input:
        vcfs=_shard_expand("results/calls/{sample}/shards/{shard}.unfiltered.vcf.gz"),
    output:
        vcf="results/calls/{sample}/{sample}.unfiltered.vcf.gz",
        idx="results/calls/{sample}/{sample}.unfiltered.vcf.gz.tbi",
    log:
        "results/logs/merge_vcfs/{sample}.log",
    params:
        args=lambda wc, input: " -I ".join(input.vcfs),
    resources:
        mem_mb=6000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx5g" MergeVcfs -I {params.args} -O {output.vcf}
        """


rule merge_stats:
    input:
        stats=_shard_expand("results/calls/{sample}/shards/{shard}.unfiltered.vcf.gz.stats"),
    output:
        "results/calls/{sample}/{sample}.merged.stats",
    log:
        "results/logs/merge_stats/{sample}.log",
    params:
        args=lambda wc, input: " -stats ".join(input.stats),
    resources:
        mem_mb=4000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx3g" MergeMutectStats -stats {params.args} -O {output}
        """


rule learn_read_orientation_model:
    # tumor cohort only -- broad's PoN WDL never sets run_orientation_bias_mixture_model_filter
    # for normals, so this step is skipped for the pon group entirely.
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        f1r2=_shard_expand("results/calls/{sample}/shards/{shard}.f1r2.tar.gz"),
    output:
        "results/calls/{sample}/{sample}.artifact-priors.tar.gz",
    log:
        "results/logs/learn_read_orientation_model/{sample}.log",
    params:
        args=lambda wc, input: " -I ".join(input.f1r2),
    resources:
        mem_mb=7000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx6g" LearnReadOrientationModel -I {params.args} -O {output}
        """


rule merge_pileup_summaries:
    # tumor cohort only -- see get_pileup_summaries.
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        tables=_shard_expand("results/calls/{sample}/shards/{shard}.pileups.table"),
        dict_=REF_DICT,
    output:
        "results/calls/{sample}/{sample}.pileups.merged.tsv",
    log:
        "results/logs/merge_pileup_summaries/{sample}.log",
    params:
        args=lambda wc, input: " -I ".join(input.tables),
    resources:
        mem_mb=4000,
    shell:
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx3g" GatherPileupSummaries \
            --sequence-dictionary {input.dict_} \
            -I {params.args} \
            -O {output}
        """


rule calculate_contamination:
    # tumor cohort only -- see get_pileup_summaries.
    wildcard_constraints:
        sample=TUMOR_SAMPLE_PATTERN,
    input:
        tumor_pileups="results/calls/{sample}/{sample}.pileups.merged.tsv",
    output:
        contamination="results/calls/{sample}/{sample}.contamination.table",
        segments="results/calls/{sample}/{sample}.segments.table",
    log:
        "results/logs/calculate_contamination/{sample}.log",
    resources:
        mem_mb=4000,
    shell:
        # Tumor-only: no -matched normal pileups.
        r"""
        exec > {log} 2>&1
        gatk --java-options "-Xmx3g" CalculateContamination \
            -I {input.tumor_pileups} \
            -O {output.contamination} \
            --tumor-segmentation {output.segments}
        """


def contamination_inputs(wildcards):
    # tumor only -- pon normals never get contamination-estimated (see module docstring)
    if is_tumor(wildcards):
        return {
            "contamination": f"results/calls/{wildcards.sample}/{wildcards.sample}.contamination.table",
            "segments": f"results/calls/{wildcards.sample}/{wildcards.sample}.segments.table",
        }
    return {}


def priors_input(wildcards):
    # tumor only -- pon normals never run LearnReadOrientationModel (see module docstring)
    if is_tumor(wildcards):
        return f"results/calls/{wildcards.sample}/{wildcards.sample}.artifact-priors.tar.gz"
    return []


def filter_extra_args(wildcards, input):
    args = []
    if is_tumor(wildcards):
        args.append(f"--contamination-table {input.contamination}")
        args.append(f"--tumor-segmentation {input.segments}")
        args.append(f"--ob-priors {input.priors}")
    extra = TUMOR_CFG.get("m2_extra_filtering_args", "") or ""
    if extra:
        args.append(extra)
    return " ".join(args)


rule filter_mutect_calls:
    input:
        unpack(lambda wc: {
            "vcf": f"results/calls/{wc.sample}/{wc.sample}.unfiltered.vcf.gz",
            "ref": REF_FASTA, "fai": REF_FAI, "dict_": REF_DICT,
            "stats": f"results/calls/{wc.sample}/{wc.sample}.merged.stats",
            "priors": priors_input(wc),
            **contamination_inputs(wc),
        }),
    output:
        vcf="results/calls/{sample}/{sample}.filtered.vcf.gz",
        idx="results/calls/{sample}/{sample}.filtered.vcf.gz.tbi",
        stats="results/calls/{sample}/{sample}.filtering.stats",
    log:
        "results/logs/filter_mutect_calls/{sample}.log",
    params:
        extra_args=filter_extra_args,
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
            --filtering-stats {output.stats} {params.extra_args}
        """


# ---------------------------------------------------------------------------
# Panel of Normals construction (skipped entirely if pon.existing_pon_vcf is set)
# ---------------------------------------------------------------------------

if BUILD_PON:

    checkpoint split_intervals_panel:
        input:
            ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT, targets=TARGETS_BED,
        output:
            directory("results/pon/panel_intervals"),
        log:
            "results/logs/split_intervals_panel/split_intervals_panel.log",
        params:
            scatter=PON_CFG["panel_scatter_count"],
            min_contig=PON_CFG["min_contig_size"],
            masked=MASKED_ARG,
        resources:
            mem_mb=6000,
        shell:
            r"""
            exec > {log} 2>&1
            mkdir -p {output}
            gatk --java-options "-Xmx5g" SplitIntervals \
                -R {input.ref} \
                -L {input.targets} \
                {params.masked} \
                --dont-mix-contigs \
                --min-contig-size {params.min_contig} \
                -scatter {params.scatter} \
                -O {output}
            """

    def panel_shard_ids(_wildcards=None):
        ck = checkpoints.split_intervals_panel.get()
        (shards,) = glob_wildcards(os.path.join(ck.output[0], "{shard}.interval_list"))
        return sorted(shards)

    rule genomicsdb_and_panel:
        input:
            vcfs=expand("results/calls/{sample}/{sample}.filtered.vcf.gz", sample=list(NORMAL_SAMPLES)),
            tbis=expand("results/calls/{sample}/{sample}.filtered.vcf.gz.tbi", sample=list(NORMAL_SAMPLES)),
            ref=REF_FASTA, fai=REF_FAI, dict_=REF_DICT,
            gnomad=GNOMAD,
            intervals="results/pon/panel_intervals/{shard}.interval_list",
        output:
            vcf="results/pon/shards/{shard}.panel.vcf",
            idx="results/pon/shards/{shard}.panel.vcf.idx",
        log:
            "results/logs/genomicsdb_and_panel/{shard}.log",
        params:
            vcf_args=lambda wc, input: " -V ".join(input.vcfs),
        resources:
            mem_mb=9000,
        shell:
            r"""
            exec > {log} 2>&1
            # GenomicsDB writes one dir of ~100 tiny files per interval (~1M files across shards),
            # which hit "No space left on device" on the shared /project GPFS. Build it on the
            # compute node's local disk instead; only the panel VCF is written to /project.
            if [ -w /scratch ]; then SCR=/scratch; else SCR=/tmp; fi
            mkdir -p $SCR/$(id -un)
            WORK=$(mktemp -d $SCR/$(id -un)/gdb_{wildcards.shard}.XXXXXX)
            trap 'rm -rf "$WORK"' EXIT
            mkdir -p $WORK/tmp
            hostname; df -h $SCR /tmp

            gatk --java-options "-Xmx8g -Djava.io.tmpdir=$WORK/tmp" GenomicsDBImport \
                --genomicsdb-workspace-path $WORK/db \
                --tmp-dir $WORK/tmp \
                -R {input.ref} \
                -V {params.vcf_args} \
                -L {input.intervals}

            gatk --java-options "-Xmx8g -Djava.io.tmpdir=$WORK/tmp" CreateSomaticPanelOfNormals \
                -R {input.ref} \
                --germline-resource {input.gnomad} \
                -V gendb://$WORK/db \
                -O {output.vcf}
            """

    def _panel_shard_expand(wildcards):
        return expand("results/pon/shards/{shard}.panel.vcf", shard=panel_shard_ids())

    rule merge_pon:
        input:
            vcfs=_panel_shard_expand,
        output:
            vcf=PON_VCF,
            idx=PON_VCF + ".tbi",
        log:
            "results/logs/merge_pon/merge_pon.log",
        params:
            args=lambda wc, input: " -I ".join(input.vcfs),
        resources:
            mem_mb=6000,
        shell:
            r"""
            exec > {log} 2>&1
            gatk --java-options "-Xmx5g" MergeVcfs -I {params.args} -O {output.vcf}
            """
