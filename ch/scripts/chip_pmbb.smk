import csv
import os

os.makedirs("cram",exist_ok=True)
os.makedirs("pileup",exist_ok=True)
os.makedirs("mutect2_out",exist_ok=True)

REF_FA=config['reference']['fasta']
IDS_CSV=config.get('chip_dx',{}).get('ids_csv','ch_psm_matched_f3_case_control_ids.csv')
PMBB_CRAM_MAP=config.get('chip_dx',{}).get('pmbb_cram_map','pmbb_to_cram2.tsv')
INTERVALS=config.get('chip_dx',{}).get('intervals',config.get('resources',{}).get('targets_bed','xgen_plus_spikein.Covered.slop10.hg38.intervals'))
PON_VCF=config.get('chip_dx',{}).get('pon_vcf','normals_merged_for_pon.vcf.gz')
TUMOR_LOD_EMIT=config.get('chip_dx',{}).get('tumor_lod_to_emit',2.0)
TUMOR_LOD_CALL=config.get('chip_dx',{}).get('initial_tumor_lod',2.0)
MIN_AF=config.get('chip_dx',{}).get('minimum_allele_fraction',0.001)
MAX_POP_AF=config.get('chip_dx',{}).get('max_population_af',0.05)
MAX_MNP=config.get('chip_dx',{}).get('max_mnp_distance',0)

def _load_pmbb_ids(path):
    with open(path,'r') as f:
        for line in f:
            s=line.strip().strip('"')
            if s:
                yield s
def _load_pmbb_to_cram(path):
    out={}
    with open(path,'r') as f:
        for line in f:
            row=line.strip().split('\t')
            if len(row)>=2 and row[0].strip() and row[1].strip():
                out[row[0].strip()]=row[1].strip()
    return out
SAMPLES=list(_load_pmbb_ids(IDS_CSV))
PMBB_TO_CRAM=_load_pmbb_to_cram(PMBB_CRAM_MAP)


rule all:
    input:
        #expand("pileup/{sample}.pileup.tsv",sample=SAMPLES),
        expand("mutect2_out/{sample}.mutect2.pon.norm.vep.vcf",sample=SAMPLES)

# download from DNA nexus
rule download_cram:
    output:
        cram=temp("cram/{cram_id}.oqfe.cram"),
        crai=temp("cram/{cram_id}.oqfe.crai")
    resources:
        download_slot=1#io_limit is a cooler name...
    shell:
        """
        dx download -f data/CRAM/{wildcards.cram_id}.oqfe.cram -o cram/{wildcards.cram_id}.oqfe.cram
        dx download -f data/CRAM/{wildcards.cram_id}.oqfe.crai -o cram/{wildcards.cram_id}.oqfe.crai
        """

# runs gatk mutect2 in tumor only w panel of normals 
rule mutect2_chip:
    input:
        cram=lambda wildcards: "cram/"+PMBB_TO_CRAM[wildcards.sample]+".oqfe.cram",
        crai=lambda wildcards: "cram/"+PMBB_TO_CRAM[wildcards.sample]+".oqfe.crai"
    output:
        vcf="mutect2_out/{sample}.mutect2.pon.vcf.gz"
    params:
        ref=REF_FA,
        intervals=INTERVALS,
        pon=PON_VCF,
        tumor_lod_emit=TUMOR_LOD_EMIT,
        tumor_lod_call=TUMOR_LOD_CALL,
        min_af=MIN_AF,
        max_pop_af=MAX_POP_AF,
        max_mnp=MAX_MNP
    threads:
        4
    shell:
        """
        gatk Mutect2 \
        -R {params.ref} \
        -L {params.intervals} \
        -I {input.cram} \
        -pon {params.pon} \
        --max-mnp-distance {params.max_mnp} \
        --tumor-lod-to-emit {params.tumor_lod_emit} \
        --initial-tumor-lod {params.tumor_lod_call} \
        --minimum-allele-fraction {params.min_af} \
        --max-population-af {params.max_pop_af} \
        -O {output.vcf}
        """

# normalizes vcf 
rule bcftools_norm:
    input:
        vcf="mutect2_out/{sample}.mutect2.pon.vcf.gz"
    output:
        vcf="mutect2_out/{sample}.mutect2.pon.norm.vcf"
    shell:
        """
        bcftools reheader -s rename.txt {input.vcf} | bcftools norm -m-both -Ov -o {output.vcf}
        """

# runs vep 
rule vep:
    input:
        vcf="mutect2_out/{sample}.mutect2.pon.norm.vcf"
    output:
        vcf="mutect2_out/{sample}.mutect2.pon.norm.vep.vcf"
    shell:
        """
        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
        singularity run --pwd "$PWD" -B "$PWD":"$PWD" -H "$PWD":"$PWD" \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i $PWD/{input.vcf} \
        -o $PWD/{output.vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --mane --everything --canonical \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta {REF_FA} \
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

# pileup region... 
rule pileup_region:
    input:
        cram=lambda wildcards: "cram/"+PMBB_TO_CRAM[wildcards.sample]+".oqfe.cram",
        crai=lambda wildcards: "cram/"+PMBB_TO_CRAM[wildcards.sample]+".oqfe.crai"
    output:
        tsv="pileup/{sample}.pileup.tsv"
    shell:
        """
        export SINGULARITY_TMPDIR=/scratch/$USER/sing_tmp
        export SINGULARITY_CACHEDIR=/scratch/$USER/sing_cache
        singularity exec \
        --bind "$(pwd):/home/bwubb" \
        --bind /scratch:/scratch \
        --pwd /home/bwubb pileup_region_latest.sif pileup_region \
        /home/bwubb/u2af1_vars.txt \
        /home/bwubb/{input.cram} \
        /home/bwubb/resources/GRCh38_full_analysis_set_plus_decoy_hla.fa > {output}
        """
