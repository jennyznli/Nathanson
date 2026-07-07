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
        # Final outputs
        expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant_STATUS.regenie", 
               PROJECT=PROJECT, CHR=CHROMOSOMES),
        expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_STATUS.regenie", 
               PROJECT=PROJECT, CHR=CHROMOSOMES),
        f"reports/{PROJECT}_regenie_report.html"

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
        fam=config['input']['exome_chr'].format(CHR=1) + ".fam"
    output:
        f"{RUN_DIR}/preprocess/{PROJECT}.exome.all_chr.step1.id"
    shell:
        """
        awk '{{print $1}}' {input.fam} > {output}
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
        bed=lambda wildcards: config['input']['exome_chr'].format(CHR=wildcards.CHR) + ".bed",
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
        (plink2 --bfile {params.input_prefix} \
               --update-sex {params.sex_info} \
               --keep {input.samples} \
               --threads {threads} \
               --geno 0.1 \
               --double-id \
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
rule preprocess_regenie:
    input:
        vep_file=config['input']['vep_file'],
        samples=f"{RUN_DIR}/preprocess/{PROJECT}.samples.txt",
        controls=config['input']['controls'],
        cases=config['input']['cases'],
        step1_covariates=config['input']['step1_covariates'],
        step2_covariates=config['input']['step2_covariates'],
    output:
        annotation=f"{RUN_DIR}/input/{PROJECT}.regenie.annotation.txt",
        set_file=f"{RUN_DIR}/input/{PROJECT}.regenie.set.txt",
        mask=f"{RUN_DIR}/input/{PROJECT}.regenie.mask.txt",
        covar_step1=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step1.txt",
        covar_step2=f"{RUN_DIR}/input/{PROJECT}.regenie.covar.step2.txt",
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt"
    params:
        output_prefix=f"{RUN_DIR}/input/{PROJECT}.regenie"
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.preprocess_regenie.log"
    shell:
        """
        (python preprocess_regenie.py \
        --vep-file {input.vep_file} \
        --samples {input.samples} \
        --cases {input.cases} \
        --controls {input.controls} \
        --step1-covariates {input.step1_covariates} \
        --step2-covariates {input.step2_covariates} \
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
        lowmem_prefix=f"{RUN_DIR}/output/tmp_rg_"
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.step1_regenie.log"
    shell:
        """
        (regenie --step 1 \
                --pgen {params.input_basename} \
                --covarFile {input.covar} \
                --phenoFile {input.pheno} \
                --extract {input.snplist} \
                --bsize 1000 \
                --threads {threads} \
                --gz \
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
        f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant_STATUS.regenie"
    params:
        input_basename=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}",
        output_basename=f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant"
    threads: 16
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.step2_single_variant.log"
    shell:
        """
        (regenie --step 2 \
                --pgen {params.input_basename} \
                --phenoFile {input.pheno} \
                --covarFile {input.covar} \
                --bt \
                --firth \
                --approx \
                --pThresh 0.999 \
                --firth-se \
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
        f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_STATUS.regenie"
    params:
        input_basename=f"{RUN_DIR}/preprocess/{PROJECT}.chr{{CHR}}",
        output_basename=f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based"
    threads: 8
    log:
        f"{RUN_DIR}/logs/{PROJECT}.chr{{CHR}}.step2_gene_based.log"
    shell:
        """
        (regenie --step 2 \
                --pgen {params.input_basename} \
                --phenoFile {input.pheno} \
                --covarFile {input.covar} \
                --bt \
                --firth \
                --approx \
                --pThresh 0.999 \
                --firth-se \
                --pred {input.pred_list} \
                --anno-file {input.annotation} \
                --set-list {input.set_file} \
                --mask-def {input.mask} \
                --build-mask 'max' \
                --write-mask-snplist \
                --check-burden-files \
                --threads {threads} \
                --af-cc \
                --bsize 200 \
                --vc-tests skat,skato \
                --out {params.output_basename}) 2>&1 | tee {log}
        """

# Generate REGENIE analysis report
rule generate_report:
    input:
        single_variant=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_single_variant_STATUS.regenie", CHR=CHROMOSOMES),
        gene_based=expand(f"{RUN_DIR}/output/{PROJECT}.chr{{CHR}}.step2_gene_based_STATUS.regenie", CHR=CHROMOSOMES),
        pheno=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        rmd="report_regenie.Rmd"
    output:
        html=f"reports/{PROJECT}_regenie_report.html",
        csv_single=f"reports/{PROJECT}_single_variant_results.csv",
        csv_gene=f"reports/{PROJECT}_gene_based_results.csv"
    params:
        project=PROJECT,
        input_dir=f"{RUN_DIR}/output",
        output_dir="reports",
        single_variant_pattern=lambda wildcards: f"{PROJECT}.chr{{CHR}}.step2_single_variant_STATUS.regenie",
        gene_based_pattern=lambda wildcards: f"{PROJECT}.chr{{CHR}}.step2_gene_based_STATUS.regenie",
        pheno_file=f"{RUN_DIR}/input/{PROJECT}.regenie.pheno.txt",
        report_title=f"SIMPLEXO ExWAS Report - {PROJECT}",
        author=config.get('report', {}).get('author', 'Analysis Report'),
        exome_sig=config.get('report', {}).get('exome_sig_threshold', 2.5e-6),
        fdr_threshold=config.get('report', {}).get('fdr_threshold', 0.05)
    threads: 1
    log:
        f"{RUN_DIR}/logs/{PROJECT}.generate_report.log"
    shell:
        """
        Rscript -e "rmarkdown::render(
            input = '{input.rmd}',
            output_file = basename('{output.html}'),
            output_dir = '{params.output_dir}',
            params = list(
                project = '{params.project}',
                input_dir = '{params.input_dir}',
                output_dir = '{params.output_dir}',
                pheno_file = '{params.pheno_file}',
                single_variant_pattern = '{params.single_variant_pattern}',
                gene_based_pattern = '{params.gene_based_pattern}',
                report_title = '{params.report_title}',
                author = '{params.author}',
                exome_sig_threshold = {params.exome_sig},
                fdr_threshold = {params.fdr_threshold}
            )
        )" 2>&1 | tee {log}
        """
