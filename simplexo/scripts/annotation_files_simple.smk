#needs vep_vcf_parser2.py

# =============================================================================
# ANNOTATION FILES (SIMPLE) -- annotation_files_simple.smk
# =============================================================================
# Preliminary-run copy of annotation_files3.smk with rare-variant QC skipped
# entirely (its freeze2-vs-freeze3 logic needs generalizing to N freezes now
# that freeze4 exists -- see rare_variant_qc.smk / rare_variant_qc.py for the
# real version, to be reintroduced later). No rare_variant_qc.smk include, no
# exclusion step, no pathogenic_vus cross-check, no sample subsetting -- VEP
# starts straight off site-qc.bcf and goes straight to the final annotation
# pgen/vcf.
# =============================================================================

CHROMOSOMES_AUTOSOMAL=list(range(1,23))

wildcard_constraints:
    CHR='[0-9]+'

rule annotation_files_simple:
    input:
        expand("data/preprocess/chr{CHR}.annotation.pgen",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.pvar",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.psam",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.vcf",CHR=CHROMOSOMES_AUTOSOMAL),
        expand("data/preprocess/chr{CHR}.annotation.no_sample.vep.report.csv",CHR=CHROMOSOMES_AUTOSOMAL),
        final_variant_stats=expand("data/qc/reports/chr{CHR}.final.variant_types.txt",CHR=CHROMOSOMES_AUTOSOMAL),

# -----------------------------------------------------------------------------
# 1. No sample subsetting (VEP starts immediately). Output: no-sample VCF,
#    derived straight from site-qc.bcf -- downstream rules that need
#    genotypes read site-qc.bcf directly, no intermediate copy needed.
# -----------------------------------------------------------------------------
rule bcftools_no_sample_vcf:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf"
    output:
        vcf="data/bcftools/chr{CHR}.site-qc.no_sample.vcf"
    shell:
        "bcftools view -G -Ov -o {output.vcf} {input.bcf}"

# -----------------------------------------------------------------------------
# 2. VEP annotation -- written directly to its final preprocess path (no
#    rare-variant-QC exclusion step to stage/copy through anymore).
# -----------------------------------------------------------------------------
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

rule count_variant_types_final:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf"
    output:
        "data/qc/reports/chr{CHR}.final.variant_types.txt"
    shell:
        """
        echo "TYPE COUNT" > {output}
        echo "TOTAL $(bcftools view -H {input.bcf} | wc -l)" >> {output}
        echo "SNP $(bcftools view -v snps -H {input.bcf} | wc -l)" >> {output}
        echo "INDEL $(bcftools view -v indels -H {input.bcf} | wc -l)" >> {output}
        echo "MNP $(bcftools view -v mnps -H {input.bcf} | wc -l)" >> {output}
        echo "OTHER $(bcftools view -v other -H {input.bcf} | wc -l)" >> {output}
        """

# -----------------------------------------------------------------------------
# 3. plink2 pgen conversion -- directly from site-qc.bcf, no rare-variant-QC
#    exclusion applied (skipped for this preliminary run).
# -----------------------------------------------------------------------------
rule plink2_annotation_pgen:
    input:
        bcf="data/bcftools/chr{CHR}.site-qc.bcf",
        sex_file="data/plink/chrX.sex_update.txt"
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

# -----------------------------------------------------------------------------
# 4. Parse final annotation VEP to report CSV.
# -----------------------------------------------------------------------------
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
        python vep_vcf_parser2.py -i {input} -o {output} -m no_sample
        """
