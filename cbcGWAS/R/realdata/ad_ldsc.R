# =============================================================================
#  realdata_ad_ldsc.R
#  Apply the conditional-prevalence correction to LDSC liability-scale h^2
#  estimates from a published APOE-adjusted Alzheimer's disease GWAS.
# =============================================================================
#
#  WHAT THIS SCRIPT DOES
#  ---------------------
#  1. PREFLIGHT phase: verify packages, reference files, and column names
#     before committing to long downloads / runs
#  2. RUN phase: normalize columns, munge, run LDSC, apply kappa correction
#
#  PREREQUISITES
#  -------------
#  R packages:
#    install.packages("data.table")
#    devtools::install_github("GenomicSEM/GenomicSEM")
#
#  Data files (see realdata_README.md):
#    eur_w_ld_chr/      -- 1000G EUR LD scores, ~5GB unpacked
#    w_hm3.snplist      -- HapMap3 reference SNPs for munging, ~3MB
#    AD_sumstats.txt.gz -- APOE-adjusted AD GWAS summary statistics
#                         (Wightman et al. 2021 or Bellenguez et al. 2022)
#
#  USAGE
#  -----
#  Step 1: edit CONFIG below to point at your local files.
#  Step 2: run with phase = "PREFLIGHT" first to confirm setup:
#            Rscript R/realdata/ad_ldsc.R 2>&1 | tee realdata_ad_preflight.log
#  Step 3: change CONFIG$phase to "RUN" and execute the full pipeline:
#            Rscript R/realdata/ad_ldsc.R 2>&1 | tee realdata_ad.log
#
#  OUTPUT (RUN phase)
#  ------------------
#  realdata_ad_results.rds  -- saved data frame with corrected h^2 values
#  realdata_ad.log          -- full console output
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

source("R/realdata/helpers.R")
source("R/sims/cbcGWAS_replication.R")  # for kappa_score, lambda_marg, lambda_cond_gauss

# -----------------------------------------------------------------------------
# CONFIG -- edit these paths
# -----------------------------------------------------------------------------
CONFIG <- list(
  # "PREFLIGHT" prints checks and exits; "RUN" does the analysis
  phase = "PREFLIGHT",

  # ---- Source GWAS file (raw, gz-compressed) ----
  ad_sumstats_raw = "data/wightman2021_ad.txt.gz",

  # Choose one of OVERRIDES_WIGHTMAN_2021 / OVERRIDES_BELLENGUEZ_2022
  # from realdata_helpers.R, or build your own list().
  # Run sanity_check_sumstats(CONFIG$ad_sumstats_raw) first if unsure.
  column_overrides = OVERRIDES_WIGHTMAN_2021,

  # ---- Sample metadata ----
  # Wightman 2021 reports N_eff = 1,126,563 (combined direct + GWAX)
  # Bellenguez 2022: 487,511 European-ancestry plus proxies
  # Use N_constant only if the file lacks a per-SNP N column
  ad_N_constant   = NULL,

  # Sample case fraction (cases / total).
  # Wightman 2021: ~90,338 / 1,126,563 = 0.080
  # Bellenguez 2022: ~111,326 / 487,511 = 0.228
  ad_sample_prev  = 0.080,

  # Population prevalence of AD (lifetime). Standard value used in published
  # AD LDSC analyses (Wightman 2021, de la Fuente 2022, Liu 2025).
  ad_pop_prev     = 0.05,

  # Trait label
  trait_name      = "AD",

  # ---- LDSC reference files ----
  ld_path         = "ldsc_data/eur_w_ld_chr/",
  hm3_snplist     = "ldsc_data/w_hm3.snplist",

  # ---- Output ----
  ad_sumstats_norm = "data/ad_normalized.txt.gz",
  out_munged       = "AD",  # GenomicSEM appends .sumstats.gz
  out_results      = "realdata_ad_results.rds",

  # rho^2_APOE plug-in grid (Liu 2025 PLOS Genet range: 4-13%)
  rho2_grid       = c(0.04, 0.07, 0.10, 0.13)
)

# =============================================================================
#  PREFLIGHT phase
# =============================================================================

if (CONFIG$phase == "PREFLIGHT") {
  cat("\n>>> PREFLIGHT phase. Edit CONFIG$phase = 'RUN' to execute.\n\n")

  preflight(
    ld_path        = CONFIG$ld_path,
    hm3_snplist    = CONFIG$hm3_snplist,
    sumstats_files = list(AD = CONFIG$ad_sumstats_raw),
    check_packages = TRUE
  )

  cat("\n>>> Inspecting raw AD sumstats file:\n\n")
  if (file.exists(CONFIG$ad_sumstats_raw)) {
    sanity_check_sumstats(CONFIG$ad_sumstats_raw)
  }

  cat("\n>>> If everything above looks correct, set CONFIG$phase = 'RUN'.\n")
  quit(save = "no", status = 0)
}

# =============================================================================
#  RUN phase
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicSEM)
})

# -----------------------------------------------------------------------------
# Step 1: Normalize column names
# -----------------------------------------------------------------------------
cat("\n=== Step 1: normalizing column names ===\n")

if (!file.exists(CONFIG$ad_sumstats_norm)) {
  prepare_sumstats(
    input_file       = CONFIG$ad_sumstats_raw,
    output_file      = CONFIG$ad_sumstats_norm,
    column_overrides = CONFIG$column_overrides,
    N_constant       = CONFIG$ad_N_constant,
    apply_qc         = TRUE,
    maf_threshold    = 0.01
  )
} else {
  cat("Normalized file already exists; skipping\n")
}

# -----------------------------------------------------------------------------
# Step 2: Munge to LDSC format
# -----------------------------------------------------------------------------
cat("\n=== Step 2: munging summary statistics for LDSC ===\n")

if (!file.exists(paste0(CONFIG$out_munged, ".sumstats.gz"))) {
  munge(
    files       = CONFIG$ad_sumstats_norm,
    hm3         = CONFIG$hm3_snplist,
    trait.names = CONFIG$out_munged,
    info.filter = 0.9,
    maf.filter  = 0.01
  )
} else {
  cat("Munged file already exists; skipping\n")
}

# -----------------------------------------------------------------------------
# Step 3: Run single-trait LDSC heritability
# -----------------------------------------------------------------------------
cat("\n=== Step 3: running LDSC ===\n")

ldsc_out <- ldsc(
  traits          = paste0(CONFIG$out_munged, ".sumstats.gz"),
  sample.prev     = CONFIG$ad_sample_prev,
  population.prev = CONFIG$ad_pop_prev,
  ld              = CONFIG$ld_path,
  wld             = CONFIG$ld_path,
  trait.names     = CONFIG$trait_name,
  ldsc.log        = "AD_ldsc"
)

h2_liab    <- ldsc_out$S[1, 1]
h2_liab_SE <- sqrt(ldsc_out$V[1, 1])

# Manual conversion as a cross-check
K  <- CONFIG$ad_pop_prev
P  <- CONFIG$ad_sample_prev
t  <- qnorm(1 - K)
conversion_factor <- (K^2 * (1 - K)^2) / (dnorm(t)^2 * P * (1 - P))
h2_obs <- h2_liab / conversion_factor

cat(sprintf("\nObserved-scale h^2:    %.4f\n", h2_obs))
cat(sprintf("Liability-scale h^2:   %.4f (SE %.4f)\n", h2_liab, h2_liab_SE))

# -----------------------------------------------------------------------------
# Step 4: Apply the conditional-prevalence correction
# -----------------------------------------------------------------------------
# Marginal scaling biases h^2_liab upward by kappa_score^2; score-based
# scaling recovers truth: h2_liab_corrected = h2_liab_marginal / kappa^2.
# Standard error scales by the same factor.

cat("\n=== Step 4: applying conditional-prevalence correction ===\n")

results <- data.frame(
  rho2_APOE         = CONFIG$rho2_grid,
  kappa_score       = NA_real_,
  h2_liab_marginal  = NA_real_,
  h2_liab_corrected = NA_real_,
  pct_change        = NA_real_,
  h2_liab_corr_SE   = NA_real_
)

for (i in seq_along(CONFIG$rho2_grid)) {
  rho2  <- CONFIG$rho2_grid[i]
  kappa <- kappa_score(CONFIG$ad_pop_prev, rho2)

  results$kappa_score[i]       <- kappa
  results$h2_liab_marginal[i]  <- h2_liab
  results$h2_liab_corrected[i] <- h2_liab / kappa^2
  results$pct_change[i]        <- 100 * (1 / kappa^2 - 1)
  results$h2_liab_corr_SE[i]   <- h2_liab_SE / kappa^2
}

cat("\n--- Results ---\n")
print(results, row.names = FALSE, digits = 4)

cat(sprintf("\nKey numbers for the manuscript:\n"))
cat(sprintf("  Marginal-scaling estimate:    h^2_liab = %.4f (SE %.4f)\n",
            h2_liab, h2_liab_SE))
cat(sprintf("  At rho^2_APOE = 0.07:         h^2_liab = %.4f, change = %.1f%%\n",
            results$h2_liab_corrected[2], results$pct_change[2]))
cat(sprintf("  Range across rho^2 plug-ins:  %.4f to %.4f (%.1f%% to %.1f%%)\n",
            min(results$h2_liab_corrected), max(results$h2_liab_corrected),
            min(results$pct_change), max(results$pct_change)))

# -----------------------------------------------------------------------------
# Save outputs
# -----------------------------------------------------------------------------
saveRDS(list(
  config     = CONFIG,
  ldsc_out   = ldsc_out,
  h2_liab    = h2_liab,
  h2_liab_SE = h2_liab_SE,
  results    = results,
  timestamp  = Sys.time()
), file = CONFIG$out_results)

cat(sprintf("\n=== saved %s ===\n", CONFIG$out_results))
cat("Send this file (or the printed key numbers above) back to fill in §4.3.\n")
