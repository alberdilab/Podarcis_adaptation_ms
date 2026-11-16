# Supplementary Table S10: homogeneity of multivariate dispersion for the dissimilarity matrices
# used in the PERMANOVAs (betadisper + permutest, 9,999 permutations, seed 123).
# Within-group Acclimation vs FMT1 rows use a PAIRED test (same individuals sampled twice): per individual, difference in
# distance to the group centroid (betadisper distances, centroids of the two time points) and a Wilcoxon signed-rank test
# (exact when possible); the `statistic` column holds 'F = x' (permutest) or 'V = x' (signed-rank), `test` names the test.
# Run from the repository root: Rscript supplementary_tables/table_S10_dispersion.R
suppressMessages({library(tidyverse); library(vegan)})
load("data/data_03122025.Rdata"); load("data/objects_03122025.Rdata"); load("data/beta_16092025.Rdata")
sm <- sample_metadata
metrics <- c(Richness = "beta_q0n", Neutral = "beta_q1n", Phylogenetic = "beta_q1p")
lab <- c(Cold_control = "CC", Cold_intervention = "CI", Warm_control = "WC",
         Podarcis_liolepis = "P. liolepis", Podarcis_muralis = "P. muralis")
relab <- function(x) { for (k in names(lab)) x <- gsub(k, lab[[k]], x, fixed = TRUE); gsub("-", "–", x) }
res <- list()
run <- function(comparison, keep, grp) {
  for (m in names(metrics)) {
    M <- as.matrix(get(metrics[[m]])$S)
    ids <- keep$Tube_code[keep$Tube_code %in% rownames(M)]
    d <- as.dist(M[ids, ids]); dat <- keep[match(labels(d), keep$Tube_code), ]
    set.seed(123)
    pt <- permutest(betadisper(d, factor(dat[[grp]])), pairwise = TRUE, permutations = 9999)
    pr <- pt$pairwise$permuted
    res[[length(res) + 1]] <<- tibble(
      comparison = comparison, metric = m, n = nrow(dat), n_pairs = NA_integer_,
      test = "betadisper permutest", F = round(pt$tab[1, "F"], 2),
      statistic = paste0("F = ", round(pt$tab[1, "F"], 2)), p = round(pt$tab[1, "Pr(>F)"], 3),
      pairwise_p = if (length(pr) > 1) relab(paste(names(pr), sprintf("%.3f", pr), collapse = "; ")) else "")
  }
}
run_paired <- function(comparison, keep) {
  for (m in names(metrics)) {
    M <- as.matrix(get(metrics[[m]])$S)
    ids <- keep$Tube_code[keep$Tube_code %in% rownames(M)]
    d <- as.dist(M[ids, ids]); dat <- keep[match(labels(d), keep$Tube_code), ]
    bd <- betadisper(d, factor(dat$time_point))
    z <- tibble(ind = dat$individual, tp = dat$time_point, z = bd$distances) %>%
      pivot_wider(names_from = tp, values_from = z, values_fn = mean) %>%
      filter(!is.na(Acclimation), !is.na(FMT1))
    wt <- suppressWarnings(wilcox.test(z$Acclimation, z$FMT1, paired = TRUE))
    res[[length(res) + 1]] <<- tibble(
      comparison = comparison, metric = m, n = nrow(dat), n_pairs = nrow(z),
      test = "paired Wilcoxon signed-rank on centroid distances", F = NA_real_,
      statistic = paste0("V = ", unname(wt$statistic)), p = round(wt$p.value, 3), pairwise_p = "")
  }
}
run("Wild: species", sm %>% filter(time_point == "Wild"), "species")
run("Acclimation: species", sm %>% filter(time_point == "Acclimation"), "species")
run("Antibiotics: species", sm %>% filter(time_point == "Antibiotics"), "species")
run("FMT1: experimental groups", sm %>% filter(time_point == "FMT1"), "type")
run("FMT2: experimental groups", sm %>% filter(time_point == "FMT2"), "type")
for (ty in c("Cold_control", "Cold_intervention", "Warm_control"))
  run_paired(paste("Acclimation vs FMT1:", lab[[ty]]), sm %>% filter(time_point %in% c("Acclimation", "FMT1"), type == ty))
out <- bind_rows(res)
write_tsv(out, "supplementary_tables/table_S10_dispersion.tsv")
print(as.data.frame(out))
