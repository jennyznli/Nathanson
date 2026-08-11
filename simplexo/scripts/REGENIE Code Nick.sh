STEP1

#!/bin/bash
#SBATCH -p med-n16-64g
#SBATCH -J REGENIE_step1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=32G
#SBATCH -t 48:00:00

######################################################## 
##  Program Name    : CouchLab_regenie_step1.sh
##  Bioinformatician: Tao Ma
##  Date Created    : Tue, 4 Feb CST 2025
##  Last Modified   : Tue, 16 Feb CST 2025
##  Function        : Run Regenie Step1 (genome-wide LOCO predictions) for one or more phenotypes
##############
##  Step 1:  A subset of variants are used to fit a whole genome regression model that
##           captures 'a good fraction' of the phenotypic variance attributable to genetic effects
##           - typically a set of several hundred thousand common markers from array
##           - 2-stage approach:
##             1) Estimate W:  ridge regression to estimate set of JM/B predictors
##                (M=#markers, B=blocksize, J=predictors per block)
##                Replace M markers with matrix W of JM/B predictors
##             2) Estimate alpha: cross-validation step to account for correlation of predictors from 1)
##                a) 2nd level of ridge for quantitative traits
##                b) logistic ridge for binary traits
##
##    Genetic Prediction:  Z = W*alpha_hat
##    - Each column of W is associated with a chromosome -
##      can construct genetic prediction leaving one chrom out by ignoring those cols.
##
######################################################## 
set -x

if [ $# != 4 ]; then
    echo 'Usage: regenie_step1.sh <PHENOTYPE_FILE> <COVARIATE_FILE> <OUTDIR> <STUDY(Must be either "PG" or "TAP")>'
    exit
fi

################# 
## variables
pheno_input="$1"
covar_input="$2"
out_dir="$3"
study=$4
#PLINK .bed/.bim/.fam files contains 279180 genetic markers that used to fit a whole genome regression model that captures a good fraction of the phenotype variance attributable to genetic effects.
if [[ ${study} == "PG" ]]; then
    bed_files=/fslustre/qhs/ext_ma_tao_mayo_edu/pipeline/Regenie/PGv2/config/PGv2_step1.gxs
elif [[ ${study} == "TAP" ]]; then
    bed_files=/fslustre/qhs/ext_ma_tao_mayo_edu/pipeline/Regenie/Tapestry/config/TAP_step1.gxs
else
    echo 'STUDY must be either "PG" or "TAP"'
    exit
fi

echo "creating output directories"
mkdir -p ${out_dir}/01.step1_WGM

REGENIE=/research/bsi/tools/biotools/regenie/4.1/bin/regenie
#REGENIE=/research/bsi/tools/biotools/regenie/4.1/miniconda3/bin/regenie
## end
###################################################

${REGENIE} \
--step 1 \
--bed ${bed_files} \
--covarFile ${covar_input} \
--phenoFile ${pheno_input} \
--bsize 1000 \
--bt \
--cv 5 \
--out ${out_dir}/01.step1_WGM/fit_step1







STEP2



#!/bin/bash
#SBATCH -p cpu-short
#SBATCH -J REGENIE_step2
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 12
#SBATCH --mem=20G
#SBATCH -t 48:00:00
#SBATCH --array=1-24

######################################################## 
##  Program Name    : CouchLab_regenie_step2.sh
##  Bioinformatician: Tao Ma
##  Date Created    : Tue, 4 Feb CST 2025
##  Last Modified   : Tue, 18 Feb CST 2025
##  Function        : Run Regenie single-variant & gene-based analysis for a specific phenotypes
######################################################## 
set -x

if [ $# != 3 ]; then
    echo "Usage: regenie_step2.sh <PHENOTYPE_FILE> <COVARIATE_FILE> <OUTDIR>"
    exit
fi

################# 
## variables
pheno_input="$1"
covar_input="$2"
out_dir="$3"

config_dir=/research/bsi/projects/PI/tertiary/Couch_Fergus_coucf/s306689.ProjectGeneration/REGENIE/config
anno_file=${config_dir}/REGENIE_format_PV_12042025.txt
mask_file=${config_dir}/REGENIE_mask_01162025.txt
set_file=${config_dir}/REGENIE_format_PV_12042025.set.txt
pgen_dir=${config_dir}/pgen

REGENIE=/research/bsi/tools/biotools/regenie/4.0/bin/regenie
## end
###################################################

chrom=$SLURM_ARRAY_TASK_ID
if [[ ${chrom} -eq 23 ]]; then
    chrom=X
elif [[ ${chrom} -eq 24 ]]; then
    chrom=Y
fi

genetic_input=${pgen_dir}/chr${chrom}
step1_input=${out_dir}/01.step1_WGM/fit_step1_pred.list

## Single-variant Association Testing
${REGENIE} \
--step 2 \
--pgen ${genetic_input} \
--phenoFile ${pheno_input} \
--covarFile ${covar_input} \
--bt \
--minMAC 1 \
--firth --approx --pThresh 0.999 \
--firth-se \
--pred ${step1_input} \
--bsize 400 \
--af-cc \
--out ${out_dir}/02.step2_WESSNP/step2_single_variant_chr${chrom}

## Gene-based Association Testing
${REGENIE} \
--step 2 \
--pgen ${genetic_input} \
--phenoFile ${pheno_input} \
--covarFile ${covar_input} \
--bt \
--firth --pThresh 0.999 \
--firth-se \
--pred ${step1_input} \
--anno-file ${anno_file} \
--set-list ${set_file} \
--mask-def ${mask_file} \
--build-mask 'max' \
--write-mask-snplist \
--check-burden-files \
--af-cc \
--minMAC 1 \
--aaf-bins 0.001,0.005,0.01,1.00 \
--bsize 200 \
--out ${out_dir}/03.step2_gene/step2_gene_based_chr${chrom}










