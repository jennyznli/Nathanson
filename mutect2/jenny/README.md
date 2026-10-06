# Tumor calling pipeline (Mutect2, tumor-only, LPC/LSF)

Snakemake reimplementation of the GATK best-practices Mutect2 calling+filtering chain
(as in `../broad/mutect2.wdl` and `../bick`'s DNAnexus port), adapted for direct
execution on Penn's LPC cluster via LSF and environment modules -- no Docker/DNAnexus,
no VEP/VarDict (add those later). Includes Panel-of-Normals construction
(`../broad/mutect2_pon.wdl`, reimplemented), interval scatter/parallelism, orientation-bias
filtering, contamination estimation, and `FilterMutectCalls`.

## What this does

For the final tumor-calling cohort (`samples/tumor_samples.tsv`):

1. `SplitIntervals` scatters the target BED into N shards (`scatter_count` in config).
2. `Mutect2` runs per-shard, tumor-only, in parallel LSF jobs, WITH `-pon`,
   `--germline-resource`, `--max-population-af` (`--native-pair-hmm-threads`,
   `--f1r2-tar-gz` for orientation-bias, tumor-tuned `--minimum-allele-fraction` /
   `--tumor-lod-to-emit` / `--initial-tumor-lod`).
3. Per-shard VCFs/stats/f1r2 are merged (`MergeVcfs`, `MergeMutectStats`,
   `LearnReadOrientationModel`).
4. `GetPileupSummaries` (per shard) + `GatherPileupSummaries` + `CalculateContamination`
   estimate cross-sample contamination.
5. `FilterMutectCalls` applies the full statistical filter (contamination, tumor
   segmentation, orientation-bias priors, LOD stats) -> `{sample}.filtered.vcf.gz`.

If `pon.existing_pon_vcf` is empty, the pipeline **also builds a PoN first**, following
broad/mutect2_pon.wdl exactly -- which is a *shorter* chain, not the same one:

1. `SplitIntervals` (shared with the tumor cohort above).
2. `Mutect2` runs per-shard on each normal in `samples/normals.tsv`, WITHOUT `-pon`
   (doesn't exist yet), WITHOUT `--germline-resource`/`--max-population-af` (would be
   circular -- see table below), `--max-mnp-distance` forced to `0`.
3. `MergeVcfs` + `MergeMutectStats` (same as tumor).
4. **`LearnReadOrientationModel`/`GetPileupSummaries`/`CalculateContamination` are skipped
   entirely for normals** -- broad's PoN WDL never wires those inputs into the per-normal
   Mutect2 sub-workflow call, so they don't run.
5. `FilterMutectCalls` still runs (`CreatePanel` consumes `Mutect2.filtered_vcf`, not the
   raw unfiltered VCF) but WITHOUT `--contamination-table`/`--tumor-segmentation`/`--ob-priors`.
6. `GenomicsDBImport` + `CreateSomaticPanelOfNormals` (scattered separately, per-contig, at
   `pon.panel_scatter_count`) build the panel from every normal's `filtered.vcf.gz`.

This happens automatically -- just leave `existing_pon_vcf` blank and the DAG handles
ordering (every tumor sample's `Mutect2` call depends on the finished PoN file).

## Setup

1. Fill in `config.yaml`:
   - `reference.fasta/fai/dict`
   - `resources.gnomad` (full germline resource, e.g. `af-only-gnomad.hg38.vcf.gz`)
   - `resources.variants_for_contamination` (common-biallelic-SNP subset, e.g.
     `small_exac_common_3.hg38.vcf.gz` -- **not** the same file as `gnomad` above)
   - `resources.targets_bed`
   - `modules.gatk` -- exact `module load` name/version on LPC
   - `lsf.queue`
   - All resource VCFs must have accompanying indices (`.tbi`) already.
2. Fill in `samples/normals.tsv` and `samples/tumor_samples.tsv` (tab-separated,
   `sample_id<TAB>path/to/reads.cram`). CRAM or BAM both work; the pipeline looks for
   `.crai`/`.bai` next to the read file automatically.
3. If you already have a PoN, set `pon.existing_pon_vcf` in config.yaml and delete/ignore
   `samples/normals.tsv` -- PoN construction is skipped entirely.

## Running

```bash
mkdir -p logs/lsf

# Always dry-run first
snakemake -n --configfile config.yaml

# Submit to LSF (Snakemake <=7.x)
snakemake --profile profiles/lsf --configfile config.yaml

# Or run locally on one node with N cores (no cluster submission)
snakemake --cores 8 --use-envmodules --configfile config.yaml
```

If your Snakemake is v8+, `--cluster`/profiles like the one here were replaced by executor
plugins -- see the comment at the top of `profiles/lsf/config.yaml` for the equivalent
`--executor cluster-generic` invocation.

## PoN vs. tumor calling: what's different

Both cohorts share the same `mutect2_scatter`/`merge_vcfs`/`merge_stats`/`filter_mutect_calls`
rule definitions (parameterized by sample group), but `get_pileup_summaries`,
`merge_pileup_summaries`, `calculate_contamination`, and `learn_read_orientation_model` are
restricted to the `tumor` cohort only (`wildcard_constraints` in `Snakefile`) -- matching
broad's `mutect2_pon.wdl` behavior exactly:

| | PoN-building normals | Final tumor samples |
|---|---|---|
| `-pon` | never (panel doesn't exist yet) | yes, once built/supplied |
| `--germline-resource` | **no** (would be circular -- `CreateSomaticPanelOfNormals` itself intersects against gnomAD when building the panel) | yes |
| `--max-population-af` | **no** (requires `--germline-resource` to mean anything) | yes, `tumor_calling.max_population_af` |
| `--max-mnp-distance` | hard-coded `0` (GenomicsDBImport can't handle MNPs) | `tumor_calling.max_mnp_distance` |
| `--minimum-allele-fraction`, `--tumor-lod-to-emit`, `--initial-tumor-lod` | same tumor-tuned values as the final cohort (sensitive detection is desirable when cataloguing artifacts for the panel too) | same |
| `LearnReadOrientationModel` | **skipped entirely** | runs, feeds `--ob-priors` |
| `GetPileupSummaries`/`CalculateContamination` | **skipped entirely** | runs, feeds `--contamination-table`/`--tumor-segmentation` |
| `FilterMutectCalls` | runs, but without `--contamination-table`/`--tumor-segmentation`/`--ob-priors` | runs with all three |

## Differences from bick's DNAnexus WDL copy (flagged for your review)

- **`masked_intervals`**: plumbed in as an optional `-XL` exclusion at both `SplitIntervals`
  steps (config: `resources.masked_intervals`, empty = no-op). Bick's copy doesn't have this;
  broad's newer upstream WDL does. You said you don't have a mask list yet -- leave it blank.
- **`FilterAlignmentArtifacts`**: intentionally **not** included (would need a
  `realignment_index_bundle` built from BWA-MEM). Can add later if you decide you want it.
- **Funcotator/annotation**: intentionally **not** included, per your instructions -- add
  VEP (as in `../brad/chip_pmbb.smk`) as a separate downstream step later.
- **`max_population_af` default (0.05)**: carried over from `../brad/chip_pmbb.smk`, more
  permissive than GATK's own Mutect2 default (0.01). Worth a second look -- it only affects
  which sites get through the *calling* stage; `FilterMutectCalls` still applies its own
  separate germline-probability filter downstream regardless of this value.
- **PoN construction**: bick's DNAnexus setup has no PoN-building workflow at all (it only
  consumes an existing PoN). This pipeline adds it, following broad's `mutect2_pon.wdl` logic.

## Notes / things to sanity-check before a big submission

- `envmodules:` directives assume `module load {modules.gatk}` puts a working `gatk` on
  `PATH`. If LPC uses a different module system quirk (e.g. requires `module load
  gatk/... samtools/...` together), edit the `envmodules:` blocks in `Snakefile` directly.
- Per-rule `mem_mb`/`threads` are reasonable starting points (matching bick/broad's WDL
  runtime defaults) but not tuned for your actual panel size/depth -- adjust in `Snakefile`
  if jobs get OOM-killed or run slower than expected.
- Run `snakemake -n` after any config change to sanity-check the DAG before submitting.
