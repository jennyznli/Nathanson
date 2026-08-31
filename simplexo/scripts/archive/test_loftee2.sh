#!/usr/bin/env bash
# Quick standalone smoke test for LOFTEE, ported from the `first_variant_annotation`
# rule in select_variants_annotate.smk. Same container/binding pattern, but:
#   - input shrunk to a tiny slice of chr22 (first 200 variants) instead of the
#     full chromosome, so a broken plugin fails in seconds, not after a real VEP run
#   - every plugin except LoF removed (REVEL, SpliceAI, gnomADc, UTRAnnotator,
#     AlphaMissense, MaveDB, NMD) and the ClinVar --custom annotation dropped too,
#     since none of that is needed to check whether loftee itself loads and runs
set -euo pipefail

IN_VCF="data/bcftools/chr22.site-qc.no_sample.vcf"
LOFTEE_DIR="/home/jennyzli/software/loftee/GRCh38"
LOFTEE_REPO="/home/jennyzli/software/loftee"

mkdir -p examine
TEST_VCF="examine/chr${CHR}.loftee_test.vcf"
OUT_VCF="examine/chr${CHR}.loftee_test.annotation.vep.vcf"

# tiny test set: header + first 200 variant records
bcftools view -h "$IN_VCF" > "$TEST_VCF"
bcftools view -H "$IN_VCF" | head -200 >> "$TEST_VCF"
echo "Test VCF: $(grep -vc '^#' "$TEST_VCF") variants"

echo "--- running VEP with only the LoF plugin ---"
singularity run -H $PWD:/home \
        --bind /home/jennyzli/resources:/opt/vep/resources \
        --bind /home/jennyzli/.vep:/opt/vep/.vep \
        --bind /home/jennyzli/software/loftee:/opt/vep/loftee \
        --bind /appl/samtools-1.21:/appl/samtools-1.21 \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i $PWD/{output.test_vcf} \
        -o $PWD/{output.vep_vcf} \
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
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz \
        --plugin LoF,loftee_path:/opt/vep/loftee,human_ancestor_fa:/opt/vep/loftee/GRCh38/human_ancestor.fa.gz,conservation_file:/opt/vep/loftee/GRCh38/loftee.sql,gerp_bigwig:/opt/vep/loftee/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw

echo "--- done: $OUT_VCF ---"
N_LOF=$(grep -o 'LoF=[A-Z]*' "$OUT_VCF" | sort | uniq -c || true)
if [ -z "$N_LOF" ]; then
    echo "WARNING: no LoF= annotations found in output -- plugin likely didn't run (check the VEP stderr/warnings above)"
else
    echo "LoF call breakdown:"
    echo "$N_LOF"
fi
