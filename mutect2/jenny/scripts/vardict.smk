import os

OUTDIR = config.get("outdir", "vardict")

os.makedirs(f"{OUTDIR}/logs", exist_ok=True)

REF_FASTA = config["reference"]["fasta"]
TARGETS_BED = config["resources"]["targets_bed"]
VARDICT_AF_THR = config.get("parameters", {}).get("vardict_af_thr", 0.001)

def _load_sample_tsv(path):
    samples = {}
    if not os.path.exists(path):
        return samples
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            samples[parts[0].strip()] = parts[1].strip()
    return samples


TUMOR_SAMPLES = _load_sample_tsv(config["samples"]["cases_tsv"])
if not TUMOR_SAMPLES:
    raise ValueError(f"No samples found in {config['samples']['cases_tsv']}")

wildcard_constraints:
    sample="|".join(__import__("re").escape(s) for s in TUMOR_SAMPLES)


rule vardict_all:
    input:
        expand(
            f"{OUTDIR}/calls/{{sample}}/{{sample}}.vcf.gz",
            sample=TUMOR_SAMPLES
        ),
        expand(
            f"{OUTDIR}/calls/{{sample}}/{{sample}}.vcf.gz.tbi",
            sample=TUMOR_SAMPLES
        )

rule run_vardict:
    input:
        bam=lambda wc: TUMOR_SAMPLES[wc.sample],
        bai=lambda wc: TUMOR_SAMPLES[wc.sample] + ".bai",
        ref=REF_FASTA,
        fai=REF_FASTA + ".fai",
        bed=TARGETS_BED
    output:
        vcf=f"{OUTDIR}/calls/{{sample}}/{{sample}}.vcf.gz",
        tbi=f"{OUTDIR}/calls/{{sample}}/{{sample}}.vcf.gz.tbi"
    threads: 4
    log:
        f"{OUTDIR}/logs/{{sample}}.run_vardict.log"
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p "$(dirname {output.vcf})"
        VARDICT_HOME="$HOME/software/VarDictJava"
        "$VARDICT_HOME/build/install/VarDict/bin/VarDict" \
            -th {threads} -G {input.ref} -f {VARDICT_AF_THR} -N {wildcards.sample} \
            -b {input.bam} -c 1 -S 2 -E 3 -g 4 {input.bed} \
        | "$VARDICT_HOME/VarDict/teststrandbias.R" \
        | "$VARDICT_HOME/VarDict/var2vcf_valid.pl" -N {wildcards.sample} -A -E -f {VARDICT_AF_THR} \
        | bgzip -c > {output.vcf}
        tabix -fp vcf {output.vcf}
        """