# =============================================================================
# SELECT VARIANTS (v3 — PMBB GoldiLocks input)
# =============================================================================
#   Per-chromosome variant QC pipeline, consuming PMBB's pre-filtered
#   GoldiLocks (GL) BCFs/VCFs. PMBB has already applied, cohort-wide
#   (70,925 samples), before this pipeline ever sees the data:
#     - split/normalized multiallelics (bcftools norm, left-aligned)
#     - genotype-level no-calling: DP<7 -> ./., GQ<20 -> ./.
#     - site-level removal: SNP sites dropped only if all het calls fail
#       Allele Balance >=15%; INDEL sites dropped only if all het calls
#       fail AB >=20%
#     - gender-discordant samples excluded from the released sample set
#   (see PMBB Release 4.0 WES documentation, Sec. 3 and Sec. 6).
#
#   This pipeline only needs to do what PMBB can't do for us because it's
#   specific to our own sample subset: subset to our cases/controls, set a
#   unique chrom_pos_ref_alt variant ID (PMBB renames *sample* IDs, not
#   variant IDs, so uniqueness here is not guaranteed upstream), recompute
#   INFO tags (AC/AF/MAF/HWE/F_MISSING/...) on that subset, apply a
#   subset-specific call-rate filter, keep PASS-only sites, and convert to
#   PLINK2 pgen/pvar/psam per autosomal chromosome. Sex for chrX ploidy
#   handling comes from config['input']['sex_info'] (no local genetic-sex
#   inference — PMBB already excludes gender-discordant samples).
#
# Configs:
#  - config/select_variants.yaml
#
# Outputs:
#   - data/bcftools/chr{CHR}.site-qc.bcf
#   - data/plink/chr{CHR}.site-qc.{pgen,pvar,psam}
# =============================================================================

import os

os.makedirs('logs/lsf',exist_ok=True)
os.makedirs('data/bcftools',exist_ok=True)
os.makedirs('data/plink',exist_ok=True)
os.makedirs('data/qc/reports',exist_ok=True)

# Extract configuration values
PROJECT_NAME=config['project']['name']
CHROMOSOMES_AUTOSOMAL=list(range(1,23))

# Individual chromosome files (REQUIRED) — already PMBB GoldiLocks (GL)
# filtered: split/normalized, genotype-level DP/GQ no-called, AB site-removed.
BCF_INPUT=config['input']['chromosomes']

wildcard_constraints:
    CHR='[0-9]+'

# Top-level, driving the whole DAG
rule select_variants:
    input:
        expand("data/qc/reports/chr{CHR}.prefilter.variant_types.select_variants.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.site-qc.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.site-qc.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.site-qc.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/bcftools/chr{CHR}.site-qc.bcf",CHR=CHROMOSOMES_AUTOSOMAL)

# Subset input BCF to the sample set, set a chrom_pos_ref_alt variant ID, and
# add subset-specific INFO tags. No bcftools norm here — PMBB's GL input is
# already split/normalized upstream (confirmed in PMBB Release 4.0 WES docs).
rule bcftools_samples_region_tags:
    input:
        bcf=lambda wildcards: BCF_INPUT.format(CHR=wildcards.CHR)
    output:
        bcf="data/bcftools/chr{CHR}.input.bcf"
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
        # write out as BCF with CSI index
        bcftools view -Ob -W=csi -o {output.bcf}
        """

######### INITIAL QC #########
# (No local sex inference: PMBB already excludes gender-discordant samples
# from the released sample set. Sex for chrX ploidy handling comes from
# config['input']['sex_info'] directly — see bcf_to_plink_per_chr below.)

# Generate stats on the tagged/subset BCF (pre site-level QC).
rule bcftools_stats_prefilter:
    input:
        bcf="data/bcftools/chr{CHR}.input.bcf"
    output:
        stats="data/qc/reports/chr{CHR}.prefilter.stats"
    shell:
        "bcftools stats -s - {input.bcf} > {output.stats}"

# Count variant types before site-level filtering
rule count_variant_types_prefilter:
    input:
        bcf="data/bcftools/chr{CHR}.input.bcf"
    output:
        "data/qc/reports/chr{CHR}.prefilter.variant_types.select_variants.txt"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# Apply site-level QC: soft-filter sites with call rate below threshold as
# LowCallRate (F_MISSING recomputed on our sample subset above — PMBB's own
# GL genotype-level/AB filtering already ran cohort-wide before we ever
# subsetted, so it can't see subset-specific missingness), then keep only
# PASS sites. No NoHQHet filter here — PMBB's GoldiLocks DP<7/GQ<20
# no-calling plus SNP/indel-specific Allele Balance site removal already
# supersedes it.
rule bcftools_filter_per_chr:
    input:
        "data/bcftools/chr{CHR}.input.bcf"
    output:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
    shell:
        """
        bcftools filter -s LowCallRate -e 'INFO/F_MISSING > 0.05' -m + {input} |
        bcftools view -f PASS -W=csi -Ob -o {output.bcf}
        """

# Generate stats after bcftools filtering
rule bcftools_stats_postfilter:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.postfilter.stats"
    shell:
        "bcftools stats -s - {input.bcf} > {output}"

# Count variant types after bcftools filtering
rule count_variant_types_postfilter:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi"
    output:
        "data/qc/reports/chr{CHR}.postfilter.variant_types.txt"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# Convert filtered BCF to PLINK format. Psam has FID=IID (double ID) and SEX;
# sex updated from config['input']['sex_info'] (headerless FID/IID/SEX file,
# 1=male/2=female — plink2's --update-sex default column layout matches this
# with no header needed).
rule bcf_to_plink_per_chr:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        csi="data/bcftools/chr{CHR}.site-qc.bcf.csi",
        sex_file=config['input']['sex_info']
    output:
        pgen="data/plink/chr{CHR}.site-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.psam"
    params:
        tmp="data/plink/chr{CHR}.site-qc.temp",
        outputname="data/plink/chr{CHR}.site-qc"
    shell:
        """
        # psam with FID=IID and SEX=0 (unknown); real sex applied below via --update-sex
        bcftools query -l {input.bcf} | awk 'BEGIN {{OFS="\\t"}} {{print $1,$1,"0"}}' | cat <(echo -e "#FID\\tIID\\tSEX") - > {params.tmp}.psam
        plink2 --bcf {input.bcf} --psam {params.tmp}.psam --vcf-half-call reference --make-pgen --out {params.tmp}
        plink2 --pfile {params.tmp} --update-sex {input.sex_file} --make-pgen --out {params.outputname}
        """
