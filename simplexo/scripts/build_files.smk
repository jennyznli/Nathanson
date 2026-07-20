# =============================================================================
# BUILD FILES: pipeline to produce the regenie step1 predictor panel + sample QC.
# Nothing here feeds regenie step2 — step2 uses the raw exome files directly.
# The only job of this pipeline is: (a) QC/exclude bad samples, (b) build a
# mildly-pruned common-SNP panel for regenie step1's null model.
# =============================================================================
# Pipeline order:
#   1. Variant QC: geno, maf, hwe, SNPs only per chr → chr*.site-qc.var-qc.* (no LD on chr)
#   2. Sample missingness on site-qc.var-qc per chr → .smiss
#   3. Aggregate missingness + apply mind → high_missingness.txt (exclude list)
#   4. Mind-passing keep list → mind_passing_keep.txt (samples to keep)
#   5. Per-chromosome, mind-passing-only: prune to TWO panels in parallel —
#      strict (r2=0.2/window 200/step 20, for het QC) and step1 (r2=0.5/window 1000/
#      step 100, for the regenie null model) — then extract each → chr{CHR}.strict.*,
#      chr{CHR}.step1.* (no long-LD-region exclusion yet; long_ld_regions.bed not available;
#      all autosomes incl. chr8 included)
#   6. Merge each panel's per-chr pruned files across chromosomes → merge_strict_pruned.*,
#      merge_step1_pruned.* (mind-passing samples, pre-het-exclusion)
#   7. Heterozygosity on merge_strict_pruned → merge_strict_pruned.het
#   8. Identify het outliers from merge_strict_pruned.het → het_exclusions.txt
#   9. passing_samples.txt = mind-passing minus het outliers (plain sample list, no genotype merge)
#  10. build_step1 = merge_step1_pruned minus het outliers → data/preprocess/build_step1.*
#      (final regenie step1 null-model input)
# PCA is NOT built here — moved to a separate script so it doesn't block this pipeline.
# Regenie covariates are built separately in R from RGC's external no1KG PCA.
# =============================================================================
#needs aggregate_missingness.py
#needs het_outliers.py
# Include from main Snakefile: include: "build_files.smk"
# Expects: config, CHROMOSOMES_AUTOSOMAL, GENO_THR, MAF_THR, HWE_THR, MIND_THR.
# Main workflow provides: data/plink/chr{CHR}.site-qc.pgen, data/plink/chrX.sex_update.txt.
# =============================================================================

CHROMOSOMES_AUTOSOMAL=list(range(1,23))
GENO_THR=config.get('qc',{}).get('geno_thr',0.01)
MAF_THR=config.get('qc',{}).get('maf_thr',0.01)
HWE_THR=config.get('qc',{}).get('hwe_thr',1e-6)
MIND_THR=config.get('qc',{}).get('mind_thr',0.05)# 0.01

wildcard_constraints:
    CHR='[0-9]+'


rule build_files:
    input:
        samples="data/qc/passing_samples.txt",
        step1_pgen="data/preprocess/build_step1.pgen",
        step1_pvar="data/preprocess/build_step1.pvar",
        step1_psam="data/preprocess/build_step1.psam",
        step1_snplist="data/preprocess/build_step1.snplist",
        prefilter_stats=expand("data/qc/reports/chr{CHR}.prefilter.variant_types.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        postfilter_stats=expand("data/qc/reports/chr{CHR}.postfilter.variant_types.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        postplink_counts=expand("data/qc/reports/chr{CHR}.postplink.variant_count.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        plink_filter_metrics=expand("data/qc/reports/chr{CHR}.plink_filter_metrics.txt", CHR=CHROMOSOMES_AUTOSOMAL),
        variant_afreq=expand("data/plink/chr{CHR}.site-qc.var-qc.afreq", CHR=CHROMOSOMES_AUTOSOMAL),
        variant_vmiss=expand("data/plink/chr{CHR}.site-qc.var-qc.vmiss", CHR=CHROMOSOMES_AUTOSOMAL),
        variant_hardy=expand("data/plink/chr{CHR}.site-qc.var-qc.hardy", CHR=CHROMOSOMES_AUTOSOMAL),

# -----------------------------------------------------------------------------
# 1. Variant QC: geno, maf, hwe, SNPs only per chr. No LD pruning on chr files.
#    Output: chr{CHR}.site-qc.var-qc.pgen/pvar/psam
# -----------------------------------------------------------------------------
rule plink2_filter_variants:
    input:
        pgen="data/plink/chr{CHR}.site-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.psam"
    output:
        pgen="data/plink/chr{CHR}.site-qc.var-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.var-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.var-qc.psam",
        log="data/plink/chr{CHR}.site-qc.var-qc.log"
    params:
        inputname="data/plink/chr{CHR}.site-qc",
        outputname="data/plink/chr{CHR}.site-qc.var-qc",
        geno=GENO_THR,
        maf=MAF_THR,
        hwe=HWE_THR
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_filter_variants.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_filter_variants.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} \
        --snps-only \
        --geno {params.geno} \
        --maf {params.maf} \
        --hwe {params.hwe} midp \
        --threads {threads} \
        --make-pgen \
        --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# QC report inputs (prefilter, postfilter, postplink, plink_filter_metrics, afreq/vmiss/hardy).
# Same paths as report expects; built from build's pvar/log only (no main pipeline).
# -----------------------------------------------------------------------------
rule count_variant_types_prefilter_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.site-qc.pvar"
    output:
        "data/qc/reports/chr{CHR}.prefilter.variant_types.txt"
    resources:
        lsf_err="logs/lsf/count_variant_types_prefilter_build.chr{CHR}.e",
        lsf_out="logs/lsf/count_variant_types_prefilter_build.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        awk -v out={output} '$0 !~ /^#/ && NF>=5 {{tot++; ref=length($4); alt=length($5); if(ref==1 && alt==1) snp++; else indel++}} END {{print "TOTAL", tot+0 >> out; print "SNP", snp+0 >> out; print "INDEL", indel+0 >> out; print "OTHER", 0 >> out}}' {input.pvar}
        """

rule count_variant_types_postfilter_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.site-qc.var-qc.pvar"
    output:
        "data/qc/reports/chr{CHR}.postfilter.variant_types.txt"
    resources:
        lsf_err="logs/lsf/count_variant_types_postfilter_build.chr{CHR}.e",
        lsf_out="logs/lsf/count_variant_types_postfilter_build.chr{CHR}.o"
    shell:
        """
        echo "TYPE COUNT" > {output}
        awk -v out={output} '$0 !~ /^#/ && NF>=5 {{tot++; ref=length($4); alt=length($5); if(ref==1 && alt==1) snp++; else indel++}} END {{print "TOTAL", tot+0 >> out; print "SNP", snp+0 >> out; print "INDEL", indel+0 >> out; print "OTHER", 0 >> out}}' {input.pvar}
        """

rule count_variants_postplink_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pvar="data/plink/chr{CHR}.site-qc.var-qc.pvar"
    output:
        "data/qc/reports/chr{CHR}.postplink.variant_count.txt"
    resources:
        lsf_err="logs/lsf/count_variants_postplink_build.chr{CHR}.e",
        lsf_out="logs/lsf/count_variants_postplink_build.chr{CHR}.o"
    shell:
        """
        TOTAL=$(awk '$0 !~ /^#/' {input.pvar} | wc -l)
        echo "CHR VARIANTS" > {output}
        echo "{wildcards.CHR} $TOTAL" >> {output}
        """

rule parse_plink_filter_metrics_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        log="data/plink/chr{CHR}.site-qc.var-qc.log"
    output:
        "data/qc/reports/chr{CHR}.plink_filter_metrics.txt"
    resources:
        lsf_err="logs/lsf/parse_plink_filter_metrics_build.chr{CHR}.e",
        lsf_out="logs/lsf/parse_plink_filter_metrics_build.chr{CHR}.o"
    shell:
        """
        echo "CHR FILTER VARIANTS_REMOVED" > {output}
        GENO=$(sed -n 's/.*\\([0-9][0-9]*\\) variant(s) removed due to missing genotype data.*/\\1/p' {input.log} 2>/dev/null | head -1); GENO=${{GENO:-0}}
        MAF=$(sed -n 's/.*\\([0-9][0-9]*\\) variant(s) removed due to allele frequency threshold.*/\\1/p' {input.log} 2>/dev/null | head -1); MAF=${{MAF:-0}}
        HWE=$(sed -n 's/.*\\([0-9][0-9]*\\) variant(s) removed due to Hardy-Weinberg.*/\\1/p' {input.log} 2>/dev/null | head -1); HWE=${{HWE:-0}}
        echo "{wildcards.CHR} GENO $GENO" >> {output}
        echo "{wildcards.CHR} MAF $MAF" >> {output}
        echo "{wildcards.CHR} HWE $HWE" >> {output}
        echo "{wildcards.CHR} PRUNE 0" >> {output}
        """

rule plink2_freq_missing_hardy_build:
    wildcard_constraints:
        CHR='[0-9]+'
    input:
        pgen="data/plink/chr{CHR}.site-qc.var-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.var-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.var-qc.psam"
    output:
        afreq="data/plink/chr{CHR}.site-qc.var-qc.afreq",
        vmiss="data/plink/chr{CHR}.site-qc.var-qc.vmiss",
        hardy="data/plink/chr{CHR}.site-qc.var-qc.hardy"
    params:
        inputname="data/plink/chr{CHR}.site-qc.var-qc",
        outputname="data/plink/chr{CHR}.site-qc.var-qc"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_freq_missing_hardy_build.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_freq_missing_hardy_build.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --freq --missing --hardy midp --threads {threads} --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# 2. Sample missingness on site-qc.var-qc per chr → chr{CHR}.site-qc.var-qc.smiss
# -----------------------------------------------------------------------------
rule plink2_sample_missing_postfilter:
    input:
        pgen="data/plink/chr{CHR}.site-qc.var-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.var-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.var-qc.psam"
    output:
        smiss="data/plink/chr{CHR}.site-qc.var-qc.smiss"
    params:
        inputname="data/plink/chr{CHR}.site-qc.var-qc",
        outputname="data/plink/chr{CHR}.site-qc.var-qc"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_sample_missing_postfilter.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_sample_missing_postfilter.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --missing --threads {threads} --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# 3. Aggregate missingness across chr1-22; apply mind threshold.
#    Output: high_missingness.txt (samples to exclude), missingness_report.txt
# -----------------------------------------------------------------------------
rule aggregate_sample_missingness:
    input:
        expand("data/plink/chr{CHR}.site-qc.var-qc.smiss", CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        aggregated="data/plink/aggregated_missingness.txt",
        exclusions="data/qc/exclusions/high_missingness.txt",
        report="data/qc/reports/missingness_report.txt"
    params:
        mind_thr=MIND_THR
    resources:
        lsf_err="logs/lsf/aggregate_sample_missingness.e",
        lsf_out="logs/lsf/aggregate_sample_missingness.o"
    shell:
        """
        python aggregate_missingness.py --mind-thr {params.mind_thr} --output-aggregated {output.aggregated} --output-exclusions {output.exclusions} --output-report {output.report} {input}
        """

# -----------------------------------------------------------------------------
# 4. Mind-passing keep list: samples not in high_missingness.txt (FID IID for --keep).
#    Output: mind_passing_keep.txt
# -----------------------------------------------------------------------------
rule mind_passing_keep_samples:
    input:
        psam="data/plink/chr1.site-qc.var-qc.psam",
        exclusions="data/qc/exclusions/high_missingness.txt"
    output:
        "data/plink/mind_passing_keep.txt"
    resources:
        lsf_err="logs/lsf/mind_passing_keep_samples.e",
        lsf_out="logs/lsf/mind_passing_keep_samples.o"
    shell:
        """
        echo "#FID IID" > {output}
        join -1 2 -2 1 -v 1 <(awk 'NR>1 {{print $1,$2}}' {input.psam} | sort -k2,2) <(sort {input.exclusions}) | sort -k2,2 >> {output}
        """

# -----------------------------------------------------------------------------
# 5. Per chromosome, mind-passing samples only: prune to two independent panels
#    in parallel — strict (r2=0.2, for het QC) and step1 (r2=0.5, for the regenie
#    null model) — then extract each directly. LD windows never cross chromosome
#    boundaries, so per-chr pruning gives an identical result to pruning after
#    merging, but runs as 22 parallel jobs instead of one serial job on the full
#    merged genome, and only needs mind_passing_keep.txt (no wait on a full merge).
#    No long-LD-region exclusion yet (long_ld_regions.bed not available); all
#    autosomes incl. chr8 included.
#    Output: chr{CHR}.strict.pgen/pvar/psam, chr{CHR}.step1.pgen/pvar/psam
# -----------------------------------------------------------------------------
rule plink2_prune_chr:
    input:
        pgen="data/plink/chr{CHR}.site-qc.var-qc.pgen",
        pvar="data/plink/chr{CHR}.site-qc.var-qc.pvar",
        psam="data/plink/chr{CHR}.site-qc.var-qc.psam",
        keep="data/plink/mind_passing_keep.txt"
    output:
        strict_pgen="data/plink/chr{CHR}.strict.pgen",
        strict_pvar="data/plink/chr{CHR}.strict.pvar",
        strict_psam="data/plink/chr{CHR}.strict.psam",
        step1_pgen="data/plink/chr{CHR}.step1.pgen",
        step1_pvar="data/plink/chr{CHR}.step1.pvar",
        step1_psam="data/plink/chr{CHR}.step1.psam"
    params:
        inputname="data/plink/chr{CHR}.site-qc.var-qc",
        strict_prefix="data/plink/chr{CHR}.strict",
        step1_prefix="data/plink/chr{CHR}.step1"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_prune_chr.chr{CHR}.e",
        lsf_out="logs/lsf/plink2_prune_chr.chr{CHR}.o"
    shell:
        """
        plink2 --pfile {params.inputname} --keep {input.keep} --indep-pairwise 200 20 0.2 --threads {threads} --out {params.strict_prefix}
        plink2 --pfile {params.inputname} --keep {input.keep} --extract {params.strict_prefix}.prune.in --make-pgen --threads {threads} --out {params.strict_prefix}
        plink2 --pfile {params.inputname} --keep {input.keep} --indep-pairwise 1000 100 0.5 --threads {threads} --out {params.step1_prefix}
        plink2 --pfile {params.inputname} --keep {input.keep} --extract {params.step1_prefix}.prune.in --make-pgen --threads {threads} --out {params.step1_prefix}
        """

# -----------------------------------------------------------------------------
# 6. Merge each panel's per-chr pruned files across chromosomes.
#    Output: merge_strict_pruned.pgen/pvar/psam, merge_step1_pruned.pgen/pvar/psam
#    (mind-passing samples, pre-het-exclusion)
# -----------------------------------------------------------------------------
rule plink2_merge_strict_pruned:
    input:
        expand("data/plink/chr{CHR}.strict.pgen", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.strict.pvar", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.strict.psam", CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        pgen="data/plink/merge_strict_pruned.pgen",
        pvar="data/plink/merge_strict_pruned.pvar",
        psam="data/plink/merge_strict_pruned.psam"
    params:
        file_list="data/plink/strict_file_list.txt",
        output_prefix="data/plink/merge_strict_pruned"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_merge_strict_pruned.e",
        lsf_out="logs/lsf/plink2_merge_strict_pruned.o"
    shell:
        """
        rm -f {params.file_list}
        for chrom in {{2..22}}; do echo "data/plink/chr${{chrom}}.strict.pgen data/plink/chr${{chrom}}.strict.pvar data/plink/chr${{chrom}}.strict.psam" >> {params.file_list}; done
        plink2 --pfile data/plink/chr1.strict --pmerge-list {params.file_list} --make-pgen --double-id --threads {threads} --out {params.output_prefix}
        rm -f {params.file_list}
        """

rule plink2_merge_step1_pruned:
    input:
        expand("data/plink/chr{CHR}.step1.pgen", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.step1.pvar", CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/plink/chr{CHR}.step1.psam", CHR=CHROMOSOMES_AUTOSOMAL)
    output:
        pgen="data/plink/merge_step1_pruned.pgen",
        pvar="data/plink/merge_step1_pruned.pvar",
        psam="data/plink/merge_step1_pruned.psam"
    params:
        file_list="data/plink/step1_file_list.txt",
        output_prefix="data/plink/merge_step1_pruned"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_merge_step1_pruned.e",
        lsf_out="logs/lsf/plink2_merge_step1_pruned.o"
    shell:
        """
        rm -f {params.file_list}
        for chrom in {{2..22}}; do echo "data/plink/chr${{chrom}}.step1.pgen data/plink/chr${{chrom}}.step1.pvar data/plink/chr${{chrom}}.step1.psam" >> {params.file_list}; done
        plink2 --pfile data/plink/chr1.step1 --pmerge-list {params.file_list} --make-pgen --double-id --threads {threads} --out {params.output_prefix}
        rm -f {params.file_list}
        """

# -----------------------------------------------------------------------------
# 7. Heterozygosity on merge_strict_pruned.
#    Output: merge_strict_pruned.het
# -----------------------------------------------------------------------------
rule plink2_het:
    input:
        pgen="data/plink/merge_strict_pruned.pgen",
        pvar="data/plink/merge_strict_pruned.pvar",
        psam="data/plink/merge_strict_pruned.psam"
    output:
        "data/plink/merge_strict_pruned.het"
    params:
        inputname="data/plink/merge_strict_pruned",
        outputname="data/plink/merge_strict_pruned"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_het.e",
        lsf_out="logs/lsf/plink2_het.o"
    shell:
        """
        plink2 --pfile {params.inputname} --het --threads {threads} --out {params.outputname}
        """

# -----------------------------------------------------------------------------
# 8. Identify het outliers from merge_strict_pruned.het (F = inbreeding; flag beyond mean ± 3*SD).
#    Output: het_exclusions.txt (samples to exclude), het_report.txt
# -----------------------------------------------------------------------------
rule plink2_het_outliers:
    input:
        "data/plink/merge_strict_pruned.het"
    output:
        exclusions="data/qc/exclusions/het_exclusions.txt",
        report="data/qc/reports/het_report.txt"
    resources:
        lsf_err="logs/lsf/plink2_het_outliers.e",
        lsf_out="logs/lsf/plink2_het_outliers.o"
    shell:
        """
        python het_outliers.py {input} --output-exclusions {output.exclusions} --output-report {output.report}
        """

# -----------------------------------------------------------------------------
# 9. passing_samples.txt = mind-passing minus het outliers. Plain sample list,
#    no genotype merge needed. Assumes het_exclusions.txt has no header (it's
#    consumed directly by plink2 --remove elsewhere without header-skipping).
#    Output: data/qc/passing_samples.txt
# -----------------------------------------------------------------------------
rule build_passing_samples:
    input:
        keep="data/plink/mind_passing_keep.txt",
        het_exclusions="data/qc/exclusions/het_exclusions.txt"
    output:
        "data/qc/passing_samples.txt"
    resources:
        lsf_err="logs/lsf/build_passing_samples.e",
        lsf_out="logs/lsf/build_passing_samples.o"
    shell:
        """
        awk 'NR==FNR {{if (FNR>1) keep[$2]=1; next}} {{del[$2]=1}} END {{for (s in keep) if (!(s in del)) print s}}' {input.keep} {input.het_exclusions} | sort > {output}
        """

# -----------------------------------------------------------------------------
# 10. build_step1 = merge_step1_pruned minus het outliers. Final regenie step1
#     null-model input.
#     Output: data/preprocess/build_step1.pgen/pvar/psam/snplist
# -----------------------------------------------------------------------------
rule plink2_create_build_step1:
    input:
        pgen="data/plink/merge_step1_pruned.pgen",
        pvar="data/plink/merge_step1_pruned.pvar",
        psam="data/plink/merge_step1_pruned.psam",
        het_exclusions="data/qc/exclusions/het_exclusions.txt"
    output:
        pgen="data/preprocess/build_step1.pgen",
        pvar="data/preprocess/build_step1.pvar",
        psam="data/preprocess/build_step1.psam",
        snplist="data/preprocess/build_step1.snplist"
    params:
        inputname="data/plink/merge_step1_pruned",
        output_prefix="data/preprocess/build_step1"
    threads: 16
    resources:
        lsf_err="logs/lsf/plink2_create_build_step1.e",
        lsf_out="logs/lsf/plink2_create_build_step1.o"
    shell:
        """
        plink2 --pfile {params.inputname} --remove {input.het_exclusions} --make-pgen --write-snplist --threads {threads} --out {params.output_prefix}
        """


