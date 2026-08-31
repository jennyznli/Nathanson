# =============================================================================
# STANDALONE LOFTEE SMOKE TEST -- full plugin stack, tiny slice
# =============================================================================
#   Ported from scripts/test_loftee2.sh into a Snakemake rule so it can be run
#   the same way as the rest of the pipeline. Same full plugin stack as
#   first_variant_annotation in select_variants_annotate.smk
#   (NMD/REVEL/SpliceAI/gnomADc/UTRAnnotator/ClinVar/AlphaMissense/MaveDB)
#   plus LoF, run against a 200-variant slice instead of the full chromosome,
#   so a broken plugin -- loftee or any of the others -- fails in seconds
#   instead of after a real fork'd VEP run.
#
#   This file is intentionally separate from select_variants_annotate.smk --
#   it is not part of that pipeline's DAG. Run it explicitly:
#     snakemake -s scripts/test_loftee_full_plugins.smk --cores 16 examine/chr22.loftee_full_test.annotation.vep.vcf
#
#   Fixes vs. test_loftee2.sh:
#     - $CHR was referenced in TEST_VCF/OUT_VCF but never set (only IN_VCF
#       hardcoded chr22), so those filenames actually resolved to
#       "examine/chr.loftee_test.vcf" -- a wildcard fixes this naturally.
#     - bind path typo: test_loftee2.sh bound
#       /home/jennyzli/loftee:/home/jennyzli/loftee, but the loftee repo is
#       actually at /home/jennyzli/software/loftee -- the real dir was never
#       bind-mounted into the container.
#     - missing --dir_plugins: per the proven-working syntax (see
#       select_variants_annotate.smk's first_variant_annotation rule and
#       loftee integration notes), a bare `--plugin LoF,...` name needs
#       `--dir_plugins /plugins,/home/jennyzli/software/loftee` alongside it
#       to resolve at all. test_loftee2.sh had neither.
#
#   Note: first_variant_annotation (in select_variants_annotate.smk) still
#   has the UNCONFIRMED full-path plugin form
#   (--plugin /home/jennyzli/software/loftee/LoF.pm,...), which produced zero
#   LoF fields in a real chr11 run. This test rule uses the proven bare-name
#   + --dir_plugins form instead, to check whether that form still works
#   alongside the full plugin stack before porting the fix back.
# =============================================================================

rule test_loftee_full_plugins:
    input:
        input_vcf="examine/chr22.loftee_test.vcf"
    output:
        vep_vcf="examine/chr22.loftee_test.annotation.vep.vcf"
    threads: 16
    resources:
        mem_mb=48000
    log:
        "examine/logs/chr22.loftee_full_test.log"
    shell:
        ""
        (
        module load samtools/1.21
        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
            --bind /home/jennyzli/resources:/opt/vep/resources \
            --bind /home/jennyzli/.vep:/opt/vep/.vep \
            --bind /home/jennyzli/software/loftee:/home/jennyzli/software/loftee \
            --bind /home/jennyzli/software/loftee/LoF.pm:/plugins/LoF.pm \
            --bind /home/jennyzli/software/loftee/ancestral.pm:/plugins/ancestral.pm \
            --bind /home/jennyzli/software/loftee/context.pm:/plugins/context.pm \
            --bind /home/jennyzli/software/loftee/de_novo_donor.pl:/plugins/de_novo_donor.pl \
            --bind /home/jennyzli/software/loftee/extended_splice.pl:/plugins/extendedomg_splice.pl \
            --bind /home/jennyzli/software/loftee/gerp_dist.pl:/plugins/gerp_dist.pl \
            --bind /home/jennyzli/software/loftee/loftee_splice_utils.pl:/plugins/loftee_splice_utils.pl \
            --bind /home/jennyzli/software/loftee/splice_site_scan.pl:/plugins/splice_site_scan.pl \
            --bind /home/jennyzli/software/loftee/svm.pl:/plugins/svm.pl \
            --bind /home/jennyzli/software/loftee/TissueExpression.pm:/plugins/TissueExpression.pm \
            --bind /home/jennyzli/software/loftee/utr_splice.pl:/plugins/utr_splice.pl \
            --bind /home/jennyzli/software/loftee/maxEntScan:/plugins/maxEntScan \
            --bind /appl/samtools-1.21/bin/samtools:/appl/samtools-1.21/bin/samtools \
            /appl/containers/vep112.sif vep \
            --dir /opt/vep/.vep \
            -i "$PWD/{input.input_vcf}" \
            -o "$PWD/{output.vep_vcf}" \
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
        ) 2>&1 | tee {log}
        """
