from datetime import datetime
import os

# Get configuration values
PROJECT = config['project']['name']
RUN_DIR = config['project'].get('run_dir', '.')  # Default to current directory if not specified
CHROMOSOMES = list(range(1, 23))  # Chromosomes 1-22
date = datetime.now().strftime("%Y%m%d")

# Create run-specific directories
os.makedirs(f"{RUN_DIR}/logs", exist_ok=True)
os.makedirs(f"{RUN_DIR}/preprocess", exist_ok=True)
os.makedirs(f"{RUN_DIR}/input", exist_ok=True)
os.makedirs(f"{RUN_DIR}/output", exist_ok=True)
os.makedirs("reports", exist_ok=True)  # Shared reports directory

# Determine step 1 input type (exome or array)
STEP1_INPUT_TYPE = config['input']['step1']['input_type']

# Optional multi-phenotype support. If config['input']['pheno_file'] is set, preprocess_regenie
# reads that pre-built multi-phenotype file directly (see preprocess_regenie.py --pheno-file)
# instead of building a single STATUS column from cases/controls. config['input']['phenotypes']
# must list the exact phenotype column names present in that file, in the order regenie should
# report them -- this list is also passed to regenie as --phenoColList so the output file names
# it writes (<out>_<phenoName>.regenie) match what's declared below. Default (no pheno_file):
# unchanged single-STATUS behavior.
PHENO_FILE = config['input'].get('pheno_file')
PHENOTYPES = config['input'].get('phenotypes', ['STATUS'])

# Brad's advanced report (report tables + advanced_report_regenie.Rmd) is single-phenotype for now.
# Uses the first entry in PHENOTYPES (STATUS by default).
REPORT_PHENO = PHENOTYPES[0]

# Define conditional inputs for step 1 based on input type
def get_step1_inputs():
    if STEP1_INPUT_TYPE == 'array':
        return {
            'pgen': f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.pgen",
            'pvar': f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.pvar", 
            'psam': f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.psam",
            'snplist': f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.snplist"
        }
    else: 
        return {
            'pgen': f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.pgen",
            'pvar': f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.pvar",
            'psam': f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.psam",
            'snplist': f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.snplist"
        }

def get_step1_dependencies():
    step1_files = get_step1_inputs()
    return [step1_files['pgen'], step1_files['snplist']]

localrules: create_vep_list, create_sample_list

wildcard_constraints:
    CHR='[1-9]|1[0-9]|2[0-2]'

# Main rule
rule run_regenie:
    input:
        # Step 1 files (conditional based on input type)
        get_step1_dependencies(),
        # Step 2 files (always use exome for step 2)
        expand(f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pgen", PROJECT=PROJECT, CHR=CHROMOSOMES),
        # Final outputs (one set per phenotype in PHENOTYPES; defaults to just STATUS)
        expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant_{{PHENO}}.regenie",
               CHR=CHROMOSOMES, PHENO=PHENOTYPES),
        expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_{{PHENO}}.regenie",
               CHR=CHROMOSOMES, PHENO=PHENOTYPES),
        expand(f"reports/{PROJECT}_{{PHENO}}_regenie_report.html", PHENO=PHENOTYPES)

# Separate target: Brad's advanced report is heavier (per-chromosome plink2 geno-counts pass)
# and not part of the default `run_regenie` target -- build explicitly with:
#   snakemake --snakefile regenie.smk run_advanced_report
rule run_advanced_report:
    input:
        f"reports/{PROJECT}_advanced_regenie_report.html"

# ========== SAMPLE LIST CREATION ==========
rule create_sample_list:
    input:
        cases=config['input']['cases'],
        controls=config['input']['controls']
    output:
        samples=f"{RUN_DIR}/preprocess/{PROJECT}.samples.txt",
        samples_plink=f"{RUN_DIR}/preprocess/{PROJECT}.samples_plink.txt"
    run:
        # Read case and control IDs
        cases = set(line.strip() for line in open(input.cases) if line.strip())
        controls = set(line.strip() for line in open(input.controls) if line.strip())
        
        # Combine and sort
        all_samples = sorted(cases | controls)
        
        print(f"Cases: {len(cases)}")
        print(f"Controls: {len(controls)}")
        print(f"Total samples: {len(all_samples)}")
        
        # Write single column format
        with open(output.samples, 'w') as f:
            for sample in all_samples:
                f.write(f'{sample}\n')
        
        # Write PLINK format (FID IID, tab-delimited)
        with open(output.samples_plink, 'w') as f:
            for sample in all_samples:
                f.write(f'{sample}\t{sample}\n')

# ========== ARRAY DATA RULES (for step 1) ==========
rule array_qc:
    input:
        bed=lambda wildcards: config['input']['step1']['array_all'] + ".bed",
        keep=f"{RUN_DIR}/preprocess/{PROJECT}.samples_plink.txt"
    output:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.pgen",
        pvar=f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.pvar",
        psam=f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.psam",
        id=f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.id",
        snplist=f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.snplist"
    params:
        input_prefix=lambda wildcards: config['input']['step1']['array_all'],
        output_prefix=f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1",
        update_ids=config['input'].get('update_ids', 'input/array_update_ids.txt')
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.array.all_chr.step1.log"
    shell:
        """
        (plink2 --bfile {params.input_prefix} \
               --keep {input.keep} \
               --geno 0.1 \
               --maf 0.01 \
               --update-ids {params.update_ids} \
               --hwe 1e-6 \
               --threads {threads} \
               --double-id \
               --write-snplist \
               --write-samples \
               --no-id-header \
               --make-pgen \
               --out {params.output_prefix}) 2>&1 | tee {log}
        """
        
# ========== EXOME DATA RULES  ==========
rule get_exome_sample_ids:
    input:
        psam=config['input']['exome_chr'].format(CHR=1) + ".psam"
    output:
        f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.id"
    shell:
        """
        awk '$1!~/^#/{{print $1}}' {input.psam} > {output}
        """

# Updated to handle conditional dependency more cleanly
def get_sample_list_inputs(wildcards):
    inputs = {
        "exome_id": f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.id",
        "config_samples": f"{RUN_DIR}/preprocess/{PROJECT}.samples.txt"
    }
    if STEP1_INPUT_TYPE == 'array':
        inputs["array_id"] = f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.id"
    return inputs

rule get_sample_list:
    input:
        unpack(get_sample_list_inputs)
    output:
        f"{RUN_DIR}/preprocess/{PROJECT}.final_samples.txt"
    run:
        # Read config samples
        config_samples = set(line.strip() for line in open(input.config_samples))
        
        # Read exome samples
        exome_samples = set(line.strip().split()[0] for line in open(input.exome_id))
        
        if STEP1_INPUT_TYPE == 'array':
            # If using array for step 1, intersect array, exome, and config samples
            array_samples = set(line.strip().split()[0] for line in open(input.array_id))
            samples = sorted(array_samples & exome_samples & config_samples)
            print(f"Intersecting array ({len(array_samples)}), exome ({len(exome_samples)}), and config ({len(config_samples)}) samples")
        else:
            # If using exome for step 1, intersect exome and config samples
            samples = sorted(exome_samples & config_samples)
            print(f"Intersecting exome ({len(exome_samples)}) and config ({len(config_samples)}) samples")
        
        print(f"Final sample count: {len(samples)}")
        
        # Write in FID IID format (both columns same, tab-delimited)
        with open(output[0], 'w') as f:
            for sample in samples:
                f.write(f'{sample}\t{sample}\n')

rule exome_chr_update_sex:
    input:
        pgen=lambda wildcards: config['input']['exome_chr'].format(CHR=wildcards.CHR) + ".pgen",
        samples=f"{RUN_DIR}/preprocess/{PROJECT}.final_samples.txt"
    output:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pgen",
        pvar=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pvar",
        psam=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.psam"
    threads: 8
    params:
        sex_info=config['input']['sex_info'],
        input_prefix=lambda wildcards: config['input']['exome_chr'].format(CHR=wildcards.CHR),
        output_prefix=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}"
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.chr_update_sex.log"
    shell:
        """
        (plink2 --pfile {params.input_prefix} \
               --update-sex {params.sex_info} \
               --keep {input.samples} \
               --threads {threads} \
               --geno 0.1 \
               --make-pgen \
               --out {params.output_prefix}) 2>&1 | tee {log}
        """
        
# Filter variants for Step 1 (always define rule, only runs if needed)
rule exome_chr_step1_filter:
    input:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pgen"
    output:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1.pgen",
        pvar=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1.pvar",
        psam=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1.psam",
        snplist=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1.snplist"
    params:
        pfile=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}",
        out=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1"
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.variant_filter.log"
    shell:
        """
        (# Calculate LD-pruned list with QC filters
        plink2 --pfile {params.pfile} \
               --double-id \
               --maf 0.01 \
               --snps-only \
               --geno 0.1 \
               --hwe 1e-6 \
               --indep-pairwise 1000 100 0.5 \
               --threads {threads} \
               --out {params.out}
        
        # Extract pruned SNPs
        plink2 --pfile {params.pfile} \
               --extract {params.out}.prune.in \
               --threads {threads} \
               --double-id \
               --write-snplist \
               --make-pgen \
               --out {params.out}) 2>&1 | tee {log}
        """

# Merge all chromosomes for Step 1 (always define rule, only runs if needed)
rule exome_merge_chr_step1:
    input:
        pgen=expand(f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1.pgen", CHR=CHROMOSOMES)
    output:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.pgen",
        pvar=f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.pvar",
        psam=f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.psam"
    params:
        out=f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1",
        file_list=f"{RUN_DIR}/preprocess/file_list.txt"
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.merge_chr.log"
    shell:
        """
        (rm -f {params.file_list}
        for chrom in {{2..22}}; do 
            echo "{RUN_DIR}/preprocess/{PROJECT}.chr${{chrom}}.step1" >> {params.file_list}
        done
        plink2 --pfile {RUN_DIR}/preprocess/{PROJECT}.chr1.step1 \
               --double-id \
               --pmerge-list {params.file_list} \
               --threads {threads} \
               --make-pgen \
               --write-samples \
               --out {params.out}
        rm -f {params.file_list}) 2>&1 | tee {log}
        """

# Merge SNP lists (always define rule, only runs if needed)
rule exome_merge_chr_snplist_step1:
    input:
        expand(f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.step1.snplist", CHR=CHROMOSOMES)
    output:
        f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.snplist"
    run:
        with open(output[0], 'w') as f:
            for i in range(1, 23):
                snplist_file = f"{RUN_DIR}/preprocess/{PROJECT}.chr{i}.step1.snplist"
                if os.path.exists(snplist_file):
                    with open(snplist_file, 'r') as f2:
                        f.write(f2.read())
                else:
                    raise FileNotFoundError(f"Expected SNPlist file not found: {snplist_file}")
                    
# ========== SHARED RULES ==========
def preprocess_regenie_inputs(wildcards):
    inputs = {
        'vep_file': config['input']['vep_file'],
        'samples': f"{RUN_DIR}/preprocess/{PROJECT}.samples.txt",
        'step1_covariates': config['input']['step1_covariates'],
        'step2_covariates': config['input']['step2_covariates'],
        'gene_consequence_exclude': config['input'].get('gene_consequence_exclude', 'mask_gene_consequence_exclude.txt'),
    }
    if PHENO_FILE:
        inputs['pheno_file'] = PHENO_FILE
    else:
        inputs['cases'] = config['input']['cases']
        inputs['controls'] = config['input']['controls']
    return inputs

def preprocess_regenie_pheno_args(wildcards, input):
    if PHENO_FILE:
        return f"--pheno-file {input.pheno_file}"
    return f"--cases {input.cases} --controls {input.controls}"

# Optional: variant IDs (e.g. potential MNP-pairs) to scrub from the M4 synonymous negative-control
# mask. Not required -- preprocess_regenie.py skips cleanly with a warning if unset/missing.
NEGATIVE_CONTROL_BLACKLIST = config['input'].get('negative_control_blacklist', '')

rule preprocess_regenie:
    input:
        unpack(preprocess_regenie_inputs)
    output:
        annotation=f"{RUN_DIR}/input/{PROJECT}.regenie.annotation.txt",
        set_file=f"{RUN_DIR}/input/{PROJECT}.regenie.set.txt",
        mask=f"{RUN_DIR}/input/{PROJECT}.regenie.mask.txt",
        covar_step1=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step1.txt",
        covar_step2=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step2.txt",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt"
    params:
        output_prefix=f"{RUN_DIR}/input/{PROJECT}.regenie",
        pheno_args=preprocess_regenie_pheno_args,
        blacklist_arg=(f"--negative-control-blacklist {NEGATIVE_CONTROL_BLACKLIST}" if NEGATIVE_CONTROL_BLACKLIST else "")
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.preprocess_regenie.log"
    shell:
        """
        (python preprocess_regenie.py \
        --vep-file {input.vep_file} \
        --samples {input.samples} \
        {params.pheno_args} \
        --step1-covariates {input.step1_covariates} \
        --step2-covariates {input.step2_covariates} \
        --gene-consequence-exclude {input.gene_consequence_exclude} \
        {params.blacklist_arg} \
        -O {params.output_prefix}) 2>&1 | tee {log}
        """

# REGENIE Step 1: Fit null model
rule run_step1_regenie:
    input:
        pgen=lambda wildcards: f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.pgen" if STEP1_INPUT_TYPE == 'array' else f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.pgen",
        covar=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step1.txt",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        snplist=lambda wildcards: f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1.snplist" if STEP1_INPUT_TYPE == 'array' else f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.snplist"
    output:
        loco=f"{RUN_DIR}/output/{PROJECT}.step1_1.loco.gz",
        pred_list=f"{RUN_DIR}/output/{PROJECT}.step1_pred.list"
    params:
        input_basename=lambda wildcards: f"{RUN_DIR}/preprocess/{PROJECT}.array.all_chr.step1" if STEP1_INPUT_TYPE == 'array' else f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1",
        output_basename=f"{RUN_DIR}/output/{PROJECT}.step1",
        lowmem_prefix=f"{RUN_DIR}/output/tmp_rg_",
        pheno_col_list=",".join(PHENOTYPES)
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.step1_regenie.log"
    shell:
        """
        (regenie --step 1 \
                --pgen {params.input_basename} \
                --covarFile {input.covar} \
                --phenoFile {input.pheno} \
                --phenoColList {params.pheno_col_list} \
                --extract {input.snplist} \
                --bsize 1000 \
                --threads {threads} \
                --gz \
                --cv 5 \
                --bt \
                --lowmem \
                --lowmem-prefix {params.lowmem_prefix} \
                --out {params.output_basename}) 2>&1 | tee {log}
        """

# REGENIE Step 2: Single variant association testing
rule run_step2_single_variant:
    input:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pgen",
        covar=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step2.txt",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        pred_list=f"{RUN_DIR}/output/{PROJECT}.step1_pred.list"
    output:
        [f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant_{p}.regenie" for p in PHENOTYPES]
    params:
        input_basename=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}",
        output_basename=f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant",
        pheno_col_list=",".join(PHENOTYPES)
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.step2_single_variant.log"
    shell:
        """
        (regenie --step 2 \
                --pgen {params.input_basename} \
                --phenoFile {input.pheno} \
                --phenoColList {params.pheno_col_list} \
                --covarFile {input.covar} \
                --bt \
                --firth \
                --approx \
                --pThresh 0.999 \
                --catCovarList batch \
                --firth-se \
                --minMAC 1 \
                --pred {input.pred_list} \
                --bsize 400 \
                --threads {threads} \
                --af-cc \
                --out {params.output_basename}) 2>&1 | tee {log}
        """

# REGENIE Step 2: Gene-based burden testing
rule run_step2_gene_based:
    input:
        pgen=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pgen",
        covar=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step2.txt",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        annotation=f"{RUN_DIR}/input/{PROJECT}.regenie.annotation.txt",
        set_file=f"{RUN_DIR}/input/{PROJECT}.regenie.set.txt",
        mask=f"{RUN_DIR}/input/{PROJECT}.regenie.mask.txt",
        pred_list=f"{RUN_DIR}/output/{PROJECT}.step1_pred.list"
    output:
        regenie=[f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_{p}.regenie" for p in PHENOTYPES],
        # --write-mask-snplist output is per-mask-definition, not per-phenotype (one file per chromosome)
        masks_snplist=f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_masks.snplist"
    params:
        input_basename=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}",
        output_basename=f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based",
        pheno_col_list=",".join(PHENOTYPES)
    threads: 8
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.step2_gene_based.log"
    shell:
        """
        (regenie --step 2 \
                --pgen {params.input_basename} \
                --phenoFile {input.pheno} \
                --phenoColList {params.pheno_col_list} \
                --covarFile {input.covar} \
                --bt \
                --firth \
                --approx \
                --pThresh 0.999 \
                --minMAC 1 \
                --catCovarList batch \
                --firth-se \
                --pred {input.pred_list} \
                --anno-file {input.annotation} \
                --set-list {input.set_file} \
                --mask-def {input.mask} \
                --build-mask 'max' \
                --write-mask-snplist \
                --check-burden-files \
                --threads {threads} \
                --aaf-bins 0.001,0.005,0.01,1.00 \
                --af-cc \
                --bsize 200 \
                --vc-tests skat,skato \
                --out {params.output_basename}) 2>&1 | tee {log}
        """

# ========== ADVANCED REPORT TABLES (Brad's report_regenie.Rmd) ==========
# mask_variant_stats.py wants a list of VEP CSVs (one per line); our pipeline only has one
# combined VEP file, so this just wraps it in a one-line list.
rule create_vep_list:
    output:
        f"{RUN_DIR}/preprocess/vep_files.list"
    params:
        vep_file=config['input']['vep_file']
    run:
        with open(output[0], 'w') as f:
            f.write(f"{params.vep_file}\n")

rule mask_variant_stats:
    input:
        vep_list=f"{RUN_DIR}/preprocess/vep_files.list",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        pgen=expand(f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}.pgen", CHR=CHROMOSOMES)
    output:
        f"{RUN_DIR}/output/{PROJECT}.mask_variant_stats.tsv"
    params:
        # lambda deactivates Snakemake's automatic wildcard expansion -- {CHR} here is a literal
        # placeholder for mask_variant_stats.py's own internal .format(CHR=chrom), not a Snakemake
        # wildcard (this rule's output has no CHR wildcard to resolve it against).
        pfile_template=lambda wildcards: f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}"
    threads: 4
    log:
        f"{RUN_DIR}/logs/{PROJECT}.mask_variant_stats.log"
    shell:
        """
        (python mask_variant_stats.py \
        --vep-list {input.vep_list} \
        --pheno {input.pheno} \
        --pfile-template {params.pfile_template} \
        --threads {threads} \
        -o {output}) 2>&1 | tee {log}
        """

rule build_regenie_report_tables:
    input:
        gene=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_{REPORT_PHENO}.regenie", CHR=CHROMOSOMES),
        snplist=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_masks.snplist", CHR=CHROMOSOMES),
        stats=f"{RUN_DIR}/output/{PROJECT}.mask_variant_stats.tsv",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt"
    output:
        top=f"{RUN_DIR}/output/{PROJECT}.report.top_genes_ADD.tsv",
        top_skat=f"{RUN_DIR}/output/{PROJECT}.report.top_genes_SKAT.tsv",
        contrib=f"{RUN_DIR}/output/{PROJECT}.report.variant_contrib.tsv",
        consequence=f"{RUN_DIR}/output/{PROJECT}.report.consequence_matrix.tsv",
        gene_carriers=f"{RUN_DIR}/output/{PROJECT}.report.gene_carriers.tsv",
        chek2_carriers=f"{RUN_DIR}/output/{PROJECT}.report.chek2_carriers.tsv"
    params:
        gene_glob=f"{RUN_DIR}/output/{PROJECT}.chr*.step2_gene_based_{REPORT_PHENO}.regenie",
        snplist_glob=f"{RUN_DIR}/output/{PROJECT}.chr*.step2_gene_based_masks.snplist",
        out_prefix=f"{RUN_DIR}/output/{PROJECT}.report",
        # lambda deactivates Snakemake's automatic wildcard expansion -- see mask_variant_stats above
        pfile_template=lambda wildcards: f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}"
    threads: 4
    log:
        f"{RUN_DIR}/logs/{PROJECT}.build_regenie_report_tables.log"
    shell:
        """
        (python build_regenie_report_tables.py \
        --gene-glob '{params.gene_glob}' \
        --snplist-glob '{params.snplist_glob}' \
        --mask-stats {input.stats} \
        --pheno {input.pheno} \
        --pfile-template {params.pfile_template} \
        --threads {threads} \
        --out-prefix {params.out_prefix}) 2>&1 | tee {log}
        """

# Generate Brad's advanced report (M1-M4 masks, top genes, variant contrib, CHEK2 section, etc.)
# Single phenotype for now (REPORT_PHENO = PHENOTYPES[0]) -- not fanned out per-phenotype yet.
rule generate_advanced_report:
    input:
        single_variant=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant_{REPORT_PHENO}.regenie", CHR=CHROMOSOMES),
        gene_based=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_{REPORT_PHENO}.regenie", CHR=CHROMOSOMES),
        report_tables=[
            f"{RUN_DIR}/output/{PROJECT}.report.top_genes_ADD.tsv",
            f"{RUN_DIR}/output/{PROJECT}.report.top_genes_SKAT.tsv",
            f"{RUN_DIR}/output/{PROJECT}.report.variant_contrib.tsv",
            f"{RUN_DIR}/output/{PROJECT}.report.consequence_matrix.tsv",
            f"{RUN_DIR}/output/{PROJECT}.report.gene_carriers.tsv",
            f"{RUN_DIR}/output/{PROJECT}.report.chek2_carriers.tsv",
        ],
        mask_stats=f"{RUN_DIR}/output/{PROJECT}.mask_variant_stats.tsv",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        rmd="advanced_report_regenie.Rmd"
    output:
        html=f"reports/{PROJECT}_advanced_regenie_report.html"
    params:
        project=PROJECT,
        data_dir=f"{RUN_DIR}/output",
        output_dir="reports",
        pheno_file=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        step2_single_variant_mid=f"step2_single_variant_{REPORT_PHENO}",
        step2_gene_based_mid=f"step2_gene_based_{REPORT_PHENO}"
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.generate_advanced_report.log"
    shell:
        """
        export LD_LIBRARY_PATH=/usr/lib64:$LD_LIBRARY_PATH
        Rscript -e "rmarkdown::render(
            input = '{input.rmd}',
            output_file = basename('{output.html}'),
            output_dir = '{params.output_dir}',
            params = list(
                project_name = '{params.project}',
                data_dir = '{params.data_dir}',
                output_dir = '{params.output_dir}',
                pheno_file = '{params.pheno_file}',
                step2_single_variant_mid = '{params.step2_single_variant_mid}',
                step2_gene_based_mid = '{params.step2_gene_based_mid}'
            )
        )" 2>&1 | tee {log}
        """
