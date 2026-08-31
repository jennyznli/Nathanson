# =============================================================================
# SELECT VARIANTS + ANNOTATE (merged, PMBB GoldiLocks input)
# =============================================================================
#   Per-chromosome variant QC + VEP annotation, consuming PMBB's pre-filtered
#   GoldiLocks (GL) VCFs directly (config['input']['chromosomes'] already
#   points at by_chromosome_GL/{prefix}.gl.chr{N}.vcf.gz). PMBB has already
#   applied, cohort-wide (70,925 samples), before this pipeline ever sees the
#   data:
#     - split/normalized multiallelics (bcftools norm, left-aligned)
#     - genotype-level no-calling: DP<7 -> ./., GQ<20 -> ./.
#     - site-level removal: SNP sites dropped only if all het calls fail
#       Allele Balance >=15%; INDEL sites dropped only if all het calls
#       fail AB >=20%
#     - gender-discordant samples excluded from the released sample set
#   (see PMBB Release 4.0 WES documentation, Sec. 3 and Sec. 6).
#
#   What's left for us to do, because it's specific to our own sample subset:
#   subset to our cases/controls, set a unique chrom_pos_ref_alt variant ID
#   (PMBB renames *sample* IDs, not variant IDs — uniqueness is not
#   guaranteed upstream), recompute INFO tags (AC/AF/MAF/HWE/F_MISSING/...)
#   on that subset, apply a subset-specific call-rate filter, keep PASS-only
#   sites, then VEP-annotate and convert to PLINK2 pgen/pvar/psam per
#   autosomal chromosome for regenie.smk. Sex for chrX ploidy handling comes
#   from config['input']['sex_info'] (headerless FID/IID/SEX, 1/2 coded) —
#   no local genetic-sex inference, since PMBB already excludes
#   gender-discordant samples.
#
#   No bcftools stats / variant-type-count checkpoints — nothing downstream
#   consumes them (checked: not read by report_regenie.Rmd or anything else
#   in the live pipeline). If you want a quick sanity count at some point,
#   `bcftools view -H file.bcf | wc -l` on the fly is enough; no need for a
#   standing rule.
#
# Dependencies:
#  - scripts/vep_vcf_parser.py
#
# Configs:
#  - config/select_variants.yaml
#
# Outputs:
#   - data/bcftools/chr{CHR}.site-qc.bcf                          (BCF QC checkpoint)
#   - data/preprocess/chr{CHR}.annotation.{pgen,pvar,psam}        (-> regenie.smk, config['input']['exome_chr'])
#   - data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv
#   - data/preprocess/all_chr.vep.report.csv                      (-> preprocess_regenie.py --vep-file)
# =============================================================================

import os

os.makedirs('logs/lsf',exist_ok=True)
os.makedirs('data/bcftools',exist_ok=True)
os.makedirs('data/preprocess',exist_ok=True)

# Extract configuration values
CHROMOSOMES_AUTOSOMAL=list(range(1,23))

# Individual chromosome files (REQUIRED) — already PMBB GoldiLocks (GL)
# filtered: split/normalized, genotype-level DP/GQ no-called, AB site-removed.
# Currently .vcf.gz (bcftools reads vcf.gz/bcf interchangeably, no conversion
# needed on the way in).
BCF_INPUT=config['input']['chromosomes']

wildcard_constraints:
    CHR='[0-9]+'

# Top-level, driving the whole DAG
rule select_variants_annotate:
    input:
        expand("data/bcftools/chr{CHR}.site-qc.bcf",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        "data/preprocess/all_chr.vep.report.csv"

# Subset input to the sample set, set a chrom_pos_ref_alt variant ID,
# recompute subset-specific INFO tags, apply the call-rate filter, and keep
# PASS-only sites — all in one streamed bcftools pipe (no intermediate BCF
# write/read between subsetting and filtering). No bcftools norm here —
# PMBB's GL input is already split/normalized upstream (confirmed in PMBB
# Release 4.0 WES docs). No NoHQHet filter either — PMBB's GoldiLocks
# DP<7/GQ<20 no-calling plus SNP/indel-specific Allele Balance site removal
# already supersedes it; LowCallRate (F_MISSING recomputed on THIS sample
# subset) is the only site-level QC left to apply here.
rule bcftools_filter_per_chr:
    input:
        bcf=lambda wildcards: BCF_INPUT.format(CHR=wildcards.CHR)
    output:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
    params:
        samples=config['input']['samples'],
        region_flag=("-R " + config['input']['targets']) if config.get('input',{}).get('targets') else ""
    shell:
        """
        # -S: restrict to sample set; -a: trim ALT alleles/INFO not seen in remaining samples; region_flag (-R) optionally restricts to target sites
        bcftools view -S {params.samples} -a {params.region_flag} {input.bcf} |
        # overwrite variant ID with CHROM_POS_REF_ALT for unambiguous downstream joins
        # (PMBB's own ID renaming is sample-level only, not guaranteed unique per-variant)
        bcftools annotate --set-id '%CHROM\\_%POS\\_%REF\\_%ALT' |
        # compute per-site INFO tags on THIS sample subset — PMBB's own cohort-wide
        # tags (if any survive from upstream) reflect all 70,925 samples, not just
        # our cases/controls, so AC/AF/MAF/HWE/F_MISSING need recomputing here
        bcftools +fill-tags - -- -t VAF,AC,AC_Het,AC_Hom,AC_Hemi,AF,AN,NS,MAF,ExcHet,F_MISSING,HWE |
        bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + |
        bcftools view -f PASS -W=csi -Ob -o {output.bcf}
        """

# Strip genotypes for a lighter/faster VEP pass — VEP annotates variants,
# not samples, and this file never has genotype data to begin with.
rule bcftools_no_sample_vcf:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf"
    output:
        vcf="data/bcftools/chr{CHR}.site-qc.no_sample.vcf"
    shell:
        "bcftools view -G -Ov -o {output.vcf} {input.bcf}"

# vep annotation (threaded)
rule first_variant_annotation:
    input:
        "data/bcftools/chr{CHR}.site-qc.no_sample.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    threads: 16
    resources:
        mem_mb=48000
    shell:
        """
        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
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
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """

# convert to plink files — the only PLINK conversion in this pipeline.
# regenie.smk reads config['input']['exome_chr'] = data/preprocess/chr{CHR}.annotation,
# so this is the one PLINK output that's actually consumed downstream.
# Sex comes straight from config['input']['sex_info'] (headerless FID/IID/SEX,
# 1/2 coded) — no local genetic-sex inference.
rule plink2_annotation_pgen:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        sex_file=config['input']['sex_info']
    output:
        pgen="data/preprocess/chr{CHR}.annotation.pgen",
        pvar="data/preprocess/chr{CHR}.annotation.pvar",
        psam="data/preprocess/chr{CHR}.annotation.psam"
    params:
        output_prefix="data/preprocess/chr{CHR}.annotation"
    threads: 16
    shell:
        """
        plink2 --threads {threads} --bcf {input.bcf} --update-sex {input.sex_file} --double-id --vcf-half-call reference --make-pgen --out {params.output_prefix}
        """

# parse final annotation VEP to report CSV.
rule parse_annotation_no_sample_vep:
    input:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf"
    output:
        "data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv"
    resources:
        mem_mb=32000,
        lsf_err="logs/lsf/parse_annotation_no_sample_vep.chr{CHR}.e",
        lsf_out="logs/lsf/parse_annotation_no_sample_vep.chr{CHR}.o"
    shell:
        """
        source $(conda info --base)/etc/profile.d/conda.sh
        conda activate vep_parser
        python scripts/vep_vcf_parser.py -i {input} -o {output} -m no_sample
        """

# Concatenate all per-chromosome VEP report CSVs into one combined file
# (header taken from chr1, then all rows minus header from every chromosome).
# This is the file to point preprocess_regenie.py's --vep-file at.
rule combine_vep_reports:
    input:
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        "data/preprocess/all_chr.vep.report.csv"
    shell:
        """
        head -1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """
