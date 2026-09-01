# SIMPLEXO4

## 1. select_variants_annotate.smk

Input is PMBB's pre-filtered GoldiLocks (GL) per-chromosome VCF

- split/normalized multiallelics

- genotype-level no-calling (DP\<7/GQ\<20 -\> `./.`)

- SNP/indel-specific Allele Balance site removal

- exclusion of gender-discordant samples

**Submission**

``` bash
bsub -N -J select_variants_annotate -eo logs/select_variants_annotate.e -oo logs/select_variants_annotate.o \
    snakemake -s scripts/select_variants_annotate.smk --configfile config/select_variants.yaml --workflow-profile \
    profiles/lsf -j 22 
```

**Config**

+---------------------+------------+-------------------------------------------------------------------+
|                     |            |                                                                   |
+=====================+============+===================================================================+
| Key                 | Required?  | Description                                                       |
+---------------------+------------+-------------------------------------------------------------------+
| `input.chromosomes` | yes        | Per-chromosome file path template (`{CHR}` placeholder)           |
+---------------------+------------+-------------------------------------------------------------------+
| `input.samples`     | yes        | Input file path w all case + control IDs, one on each line        |
+---------------------+------------+-------------------------------------------------------------------+
| `input.sex_info`    | yes        | Headerless FID/IID/SEX file (1=male, 2=female) for `--update-sex` |
+---------------------+------------+-------------------------------------------------------------------+
| `input.targets`     | no         | Optional BED/region file to restrict sites via `-R`               |
+---------------------+------------+-------------------------------------------------------------------+

- **subset to sample set, set unique chrom_pos_ref_alt variant ID, recompute subset-specific INFO tags, apply call-rate filter, keep PASS-only sites -- one streamed bcftools pipe, no intermediate BCF writewait**
  - `bcftools_filter_per_chr`
  - input: per-chr GL VCF (config)
  - output: data/bcftools/chr{CHR}.site-qc.bcf
  - no NoHQHet filter, no sex-check -- PMBB's GL input already covers both (see above); LowCallRate (`F_MISSING` recomputed on THIS sample subset, not PMBB's full cohort) is the only site-level QC left to apply ourselves
- **strip genotypes for a lighter/faster VEP pass**
  - `bcftools_no_sample_vcf`
  - data/bcftools/chr{CHR}.site-qc.bcf -\> data/bcftools/chr{CHR}.site-qc.no_sample.vcf
- **VEP annotation by running VEP in a singularity container (`vep112.sif`)**
  - `first_variant_annotation`
  - data/bcftools/chr{CHR}.site-qc.no_sample.vcf -\> data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf
- **convert site-qc BCF to PLINK w/ double ID and updated sex -- the only PLINK conversion in this pipeline; `regenie.smk`'s `input.exome_chr` reads directly from here**
  - `plink2_annotation_pgen`
  - data/bcftools/chr{CHR}.site-qc.bcf -\> data/preprocess/chr{CHR}.annotation.{pgen,pvar,psam}
- **parses VEP annotated VCF to CSV (`-m no_sample`: per-variant annotation only, no genotype/zygosity columns), expects a `vep_parser` conda env to be activatable**
  - `parse_annotation_no_sample_vep`
  - data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf -\> data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv
- **combines all per-chromosome VEP report CSVs into one file**
  - `combine_vep_reports`
  - data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv -\> data/preprocess/all_chr.vep.report.csv

## 2. regenie.smk

**Submission**

``` bash
bsub -N -J array4 -eo logs/array4.e -oo logs/array4.o \
   snakemake -s scripts/regenie.smk --workflow-profile profiles/lsf -j 22 \
  --configfile config/array4.yaml
```

**Config**

+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| Key                                                | Required?                                        | Description                                                                                                                                                                                                                                              |
+====================================================+==================================================+==========================================================================================================================================================================================================================================================+
| `project.name`                                     | yes                                              | Run/output naming prefix (`{PROJECT}` in every output path)                                                                                                                                                                                              |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `project.run_dir`                                  | no (default `.`)                                 | Base directory for `preprocess/`, `input/`, `output/`, `logs/`                                                                                                                                                                                           |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.vep_file`                                   | yes                                              | Combined VEP CSV from `select_variants_annotate.smk`'s `combine_vep_reports` -\> `preprocess_regenie.py --vep-file`                                                                                                                                      |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.exome_chr`                                  | yes                                              | Per-chromosome exome pgen/psam path template -- Step 2 genotypes, and Step 1 exome path if used                                                                                                                                                          |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.pheno_file`                                 | no                                               | Pre-built multi-phenotype file (IID + phenotype columns). If set, used instead of building a single STATUS column from cases/controls                                                                                                                    |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.phenotypes`                                 | no (default `[STATUS]`)                          | Phenotype column names in `pheno_file`, in report order -- also passed to regenie as `--phenoColList`                                                                                                                                                    |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.step1.input_type`                           | yes                                              | `"array"` or `"exome"` -- picks the Step 1 genotype prep path                                                                                                                                                                                            |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.step1.array_all`                            | yes if array                                     | PLINK bfile prefix for the array QC step                                                                                                                                                                                                                 |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.cases`, `input.controls`                    | yes                                              | Always required (even in multi-phenotype mode) -- `create_sample_list` builds the base sample set from these regardless of `pheno_file`                                                                                                                  |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.step1_covariates`, `input.step2_covariates` | yes                                              | Covariate files (any columns, e.g. PCs + batch), passed straight through                                                                                                                                                                                 |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.sex_info`                                   | yes                                              | Sex update file for Step 2 exome genotypes                                                                                                                                                                                                               |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.mask_file`                                  | no                                               | Custom regenie mask definitions; falls back to built-in `M1 pathogenic` / `M2 pathogenic,vus` if omitted                                                                                                                                                 |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.gene_consequence_exclude`                   | no (default `mask_gene_consequence_exclude.txt`) | GENE+consequence exclusion rules, dropped before pathogenic/VUS filtering.                                                                                                                                                                               |
|                                                    |                                                  |                                                                                                                                                                                                                                                          |
|                                                    |                                                  | **Note:** unlike `mask_file`/`negative_control_blacklist`, this one is wired as a hard Snakemake input -- if this file doesn't exist (default or configured path), the DAG will fail to build with a missing-input error rather than skipping gracefully |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.negative_control_blacklist`                 | no                                               | Variant IDs to scrub from the M4 synonymous negative-control mask                                                                                                                                                                                        |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| `input.update_ids`                                 | no (default `input/array_update_ids.txt`)        | ID remap file for `array_qc`, only used on the array Step 1 path                                                                                                                                                                                         |
+----------------------------------------------------+--------------------------------------------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

- **union cases + controls into one sample list**
  - `create_sample_list`
  - input: cases, controls (config)
  - output: {PROJECT}.samples.txt, {PROJECT}.samples_plink.txt
- **Step 1 genotype prep**
  - selects one path to use for step 1 of Regenie depending on config:
  - `array`:
    - QC's imputed array (geno/MAF/HWE filters, ID remap)
    - imputed array (config) -\> {PROJECT}.array.all_chr.step1.{pgen,pvar,psam,snplist}
  - `exome:`
    - `get_exome_sample_ids -> get_sample_list -> exome_chr_update_sex:` intersects exome IDs with sample list, updates sex of these files
    - `exome_chr_step1_filter (LD-prune + QC, per chr) -> exome_merge_chr_step1 + exome_merge_chr_snplist_step1 (merge across chromosomes) -> {PROJECT}.exome.all_chr.step1.{pgen,pvar,psam,snplist}:` applies LD pruning and QC, merges them into genome-wide file
- **Step 2 genotype prep: update sex, restrict to final sample set, light per-chr missingness filter**
  - `exome_chr_update_sex`
  - input: data/preprocess/{CHR}.pgen (from select_variants_annotate.smk), sex_info (config)
  - output: {PROJECT}.chr{CHR}.pgen -- used directly by both step2 rules below
- **preprocess_regenie: builds all regenie-ready inputs**
  - `preprocess_regenie`
  - calls scripts/preprocess_regenie.py to build the annotation/set/mask files as well as case/control files
  - input: all_chr.vep.report.csv, samples, cases/controls (or pre-built multi-phenotype file), step1/step2 covariate files
  - output: {PROJECT}.regenie.{annotation,set,mask,covar.step1,covar.step2,pheno}.txt
- **Step 1: fit null model**
  - `run_step1_regenie`
  - input: Step 1 pgen/snplist, covar.step1, pheno
  - output: {PROJECT}.step1_1.loco.gz, {PROJECT}.step1_pred.list
- **Step 2: single-variant association test**
  - `run_step2_single_variant`
  - input: Step 2 per-chr pgen, covar.step2, pheno, step1_pred.list
  - output: {PROJECT}.chr{CHR}.step2_single_variant\_{PHENO}.regenie
- **Step 2: gene-based burden/SKAT test**
  - `run_step2_gene_based`
  - input: Step 2 per-chr pgen, covar.step2, pheno, annotation, set, mask, step1_pred.list
  - output: {PROJECT}.chr{CHR}.step2_gene_based\_{PHENO}.regenie
- **create HTML report**
  - `create_vep_list -> mask_variant_stats -> build_regenie_report_tables -> generate_advanced_report`
  - (one per phenotype -- top genes, variant contribution, consequence matrix, gene/CHEK2 carriers)
  - input: step2 single-variant + gene-based outputs, pheno file
  - output: reports/{PROJECT}\_{PHENO}\_report.html

## Dependencies

+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
|                                  |                             |                                                                                                                                 |
+==================================+=============================+=================================================================================================================================+
|                                  | Used by                     | Purpose                                                                                                                         |
+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| `vep_vcf_parser.py`              | select_variant_annotate.smk | select_variants_annotate.smk \| Parses VEP-annotated VCF -\> report CSV (`-m no_sample`: per-variant only, no genotype columns) |
+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| `preprocess_regenie.py`          | regenie.smk                 | Builds annotation/set/mask/covariate/pheno files from the VEP CSV                                                               |
+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| `mask_variant_stats.py`          | regenie.smk                 | Per-mask, per-gene carrier/variant stats for the report                                                                         |
+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| `build_regenie_report_tables.py` | regenie.smk                 | Assembles top-genes/variant-contrib/carrier tables for the report                                                               |
+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+
| `report_regenie.Rmd`             | regenie.smk                 | Renders the final per-phenotype HTML report                                                                                     |
+----------------------------------+-----------------------------+---------------------------------------------------------------------------------------------------------------------------------+

#### preprocess_regenie.py

- takes the combined VEP CSV, sample/case/control lists (or a pre-built multi-phenotype file), and Step 1/Step 2 covariate files, and formats them for regenie's burden-test framework
- **annotation file**: one row per unique variant ID, most-pathogenic annotation kept on duplicates
  - `pathogenic`: Variant.LoF_level == 1
  - `vus`: Variant.LoF_level == 2
  - `synonymous`: negative-control category pulled from Variant.Consequence (M4), excludes anything already claimed by pathogenic/vus and anything in `--negative-control-blacklist`
- **set file**: for each gene, chr + position (from its first variant) + comma-separated variant IDs
- **mask file**: default is `M1 pathogenic` / `M2 pathogenic,vus`. Can be overridden with `--mask-file` (custom mask lines `Mx category1,category2,...`; categories are validated against {pathogenic, vus, synonymous} and checked for zero-variant categories)
- **covariate files**: passed straight through from `--step1-covariates` / `--step2-covariates` (whatever columns they contain -- eigenvec, batch, etc.), filtered to the final sample set. If omitted, writes a minimal FID/IID-only file.
- **phenotype file**: single `STATUS` column (0/1) built from `--cases`/`--controls`, OR passed through from `--pheno-file` for multi-phenotype runs (`config['input']['phenotypes']` lists the column names, in the order regenie should report them)
- optional: `--gene-consequence-exclude` (drop specific GENE+consequence combos before pathogenic/VUS filtering)

## 

finish this....
