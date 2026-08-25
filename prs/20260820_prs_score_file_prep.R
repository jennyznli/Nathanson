# ============================================================
# BUILD PRS SCORE FILES
# ============================================================
# Produces the 4-column (MarkerName, EffectAllele, AltAllele, Effect) score
# files tgct/scripts/scripts/prs.sh expects. AltAllele (renamed here from
# OtherAllele/A2) is required so prs.sh's standardize_ids() can validate
# that each genotype source actually has the expected biallelic pair at a
# given position, instead of just matching on chr:pos alone.
library(here)
source("R/config.R")

# ========================
# 2021 GWAS (meta-analysis hits)
# ========================
gwas2021_raw <- read_excel(here("tgct", "ss", "meta-hits-062526.xlsx"), sheet = "PreviousHits")

# replicated hits only; Z = effect / SE gives a second scoring option
# alongside the raw effect size (beta)
gwas2021_replicated <- gwas2021_raw |>
    filter(Replicated == "Replicated") |>
    select(MarkerName, EffectAllele, OtherAllele, Effect, StdErr) |>
    mutate(
        Z = Effect / StdErr
    )

# effect (beta) scored version
gwas2021_score_effect <- gwas2021_replicated |>
    select(MarkerName, EffectAllele, OtherAllele, Effect) |>
    rename(AltAllele = OtherAllele)
write_tsv(gwas2021_score_effect, here("tgct", "data", "scorefiles", "gwas2021_hits_effect.txt"))

# Z-scored version
gwas2021_score_z <- gwas2021_replicated |>
    select(MarkerName, EffectAllele, OtherAllele, Z) |>
    rename(AltAllele = OtherAllele)
write_tsv(gwas2021_score_z, here("tgct", "data", "scorefiles", "gwas2021_hits_z.txt"))


# ========================
# 2026 GWAS (emilio's association results)
# ========================
gwas2026_raw <- read.table(here("tgct", "data", "emilio-scores-old.txt"), header = TRUE)

# effect (beta) scored version
gwas2026_score_effect <- gwas2026_raw |> select(SNP_b38, A1, A2, BETA) |>
    rename(MarkerName = SNP_b38) |>
    rename(EffectAllele = A1) |>
    rename(AltAllele = A2) |>
    rename(Effect = BETA)
write_tsv(gwas2026_score_effect, here("tgct", "data", "scorefiles", "gwas2026_hits_effect.txt"))

# Z-scored version
gwas2026_score_z <- gwas2026_raw |> select(SNP_b38, A1, A2, Z) |>
    rename(MarkerName = SNP_b38) |>
    rename(EffectAllele = A1) |>
    rename(AltAllele = A2) |>
    rename(Effect = Z)
write_tsv(gwas2026_score_z, here("tgct", "data", "scorefiles", "gwas2026_hits_z.txt"))
