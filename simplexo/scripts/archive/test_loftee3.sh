#!/usr/bin/env bash
# Standalone smoke test for LOFTEE, ported from the `first_variant_annotation`
# rule in select_variants_annotate.smk. Same container/binding pattern and now
# the SAME FULL PLUGIN STACK (NMD/REVEL/SpliceAI/gnomADc/UTRAnnotator/ClinVar/
# AlphaMissense/MaveDB/LoF) as that rule, matching test_loftee_full_plugins.smk
# -- only difference from the real pipeline rule is the input, shrunk to a
# tiny slice of chr22 (first 200 variants) instead of the full chromosome, so
# a broken plugin fails in seconds, not after a real VEP run.
#
# NOTE: mirrors the bind/mount pattern currently in test_loftee_full_plugins.smk
# (-H $PWD:/home, /home/jennyzli/resources + /home/jennyzli/.vep instead of
# bwubb's, loftee mounted at /opt/vep/loftee, no --dir_plugins). That pattern
# is UNCONFIRMED -- it differs from the empirically-proven form documented for
# first_variant_annotation (--pwd/-B/-H all pointed at $PWD, --dir_plugins
# /plugins,<loftee-dir> alongside a bare --plugin LoF,... name). See the
# flags raised alongside this file for specifics.
set -euo pipefail

IN_VCF="data/bcftools/chr22.site-qc.no_sample.vcf"
LOFTEE_DIR="/opt/vep/loftee/GRCh38"
LOFTEE_REPO="/opt/vep/loftee"

mkdir -p examine
TEST_VCF="examine/chr22.loftee_test.vcf"
OUT_VCF="examine/chr22.loftee_test.annotation.vep.vcf"

[ -f "$IN_VCF" ] || { echo "ERROR: missing $IN_VCF" >&2; exit 1; }

# tiny test set: header + first 200 variant records
bcftools view -h "$IN_VCF" > "$TEST_VCF"
bcftools view -H "$IN_VCF" | head -200 >> "$TEST_VCF"
echo "Test VCF: $(grep -vc '^#' "$TEST_VCF") variants"

echo "--- running VEP with the full plugin stack (NMD/REVEL/SpliceAI/gnomADc/UTRAnnotator/ClinVar/AlphaMissense/MaveDB/LoF) ---"
export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache

singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
    --bind /home/jennyzli/resources:/opt/vep/resources \
    --bind /home/jennyzli/.vep:/opt/vep/.vep \
    --bind /home/jennyzli/software/loftee:/home/jennyzli/software/loftee \
    --bind /home/jennyzli/software/loftee/LoF.pm:/plugins/LoF.pm \
    --bind /home/jennyzli/software/loftee/ancestral.pm:/plugins/ancestral.pm \
    --bind /home/jennyzli/software/loftee/context.pm:/plugins/context.pm \
    --bind /home/jennyzli/software/loftee/de_novo_donor.pl:/plugins/de_novo_donor.pl \
    --bind /home/jennyzli/software/loftee/extended_splice.pl:/plugins/extended_splice.pl \
    --bind /home/jennyzli/software/loftee/gerp_dist.pl:/plugins/gerp_dist.pl \
    --bind /home/jennyzli/software/loftee/loftee_splice_utils.pl:/plugins/loftee_splice_utils.pl \
    --bind /home/jennyzli/software/loftee/splice_site_scan.pl:/plugins/splice_site_scan.pl \
    --bind /home/jennyzli/software/loftee/svm.pl:/plugins/svm.pl \
    --bind /home/jennyzli/software/loftee/TissueExpression.pm:/plugins/TissueExpression.pm \
    --bind /home/jennyzli/software/loftee/utr_splice.pl:/plugins/utr_splice.pl \
    --bind /home/jennyzli/software/loftee/maxEntScan:/plugins/maxEntScan \
    --bind /home/jennyzli/bin/samtools:/appl/samtools-1.21/bin/samtools \
    /appl/containers/vep112.sif vep \
    --dir /opt/vep/.vep \
    -i "$PWD/$TEST_VCF" \
    -o "$PWD/$OUT_VCF" \
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
    --plugin LoF,loftee_path:/home/jennyzli/software/loftee,human_ancestor_fa:/home/jennyzli/software/loftee/GRCh38/human_ancestor.fa.gz,conservation_file:/home/jennyzli/software/loftee/GRCh38/loftee.sql,gerp_bigwig:/home/jennyzli/software/loftee/GRCh38/gerp_conservation_scores.homo_sapiens.GRCh38.bw \
    --cache_version 116


    --dir_plugins /plugins \
    --dir_plugins /opt/vep/loftee \


echo "--- done: $OUT_VCF ---"
N_LOF=$(grep -o 'LoF=[A-Z]*' "$OUT_VCF" | sort | uniq -c || true)
if [ -z "$N_LOF" ]; then
    echo "WARNING: no LoF= annotations found in output -- plugin likely didn't run (check the VEP stderr/warnings above)"
else
    echo "LoF call breakdown:"
    echo "$N_LOF"
fi
