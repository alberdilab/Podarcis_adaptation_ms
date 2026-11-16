# Are MAGs associated with the individual they were assembled from?
# Individual-level co-assemblies pool wild, acclimation, antibiotic and FMT samples
# of one lizard. This script tests (1) whether bins behave as single populations across
# those time points (contig-level coverage homogeneity, calibrated with artificial
# chimeras), (2) whether MAGs are preferentially detected in their source individual
# and (3) whether samples cluster by individual across time points.
# Run from the repository root: Rscript mag_individual_association/mag_individual_association.R

suppressMessages({library(tidyverse); library(vegan)})
set.seed(20261115)
out <- "mag_individual_association"
load("data/data_03122025.Rdata")

time_levels <- c("Wild","Acclimation","Antibiotics","Transplant1","Transplant2","FMT1","FMT2")
smeta <- sample_metadata %>% mutate(time_point = factor(time_point, levels = time_levels))

# ---- 1. contig-level coverage homogeneity ------------------------------------
cc <- read_tsv("data/contig.count.tsv.xz", show_col_types = FALSE)
contig_id <- cc$sequence_id
mat <- as.matrix(cc[,-1])
samp <- sub("_.*", "", colnames(mat))               # sum libraries of the same sample
mat <- sapply(split(seq_along(samp), samp), function(i) rowSums(mat[, i, drop = FALSE]))
rownames(mat) <- contig_id
mat <- mat[, colnames(mat) %in% smeta$Tube_code, drop = FALSE]
mag_of <- sub("@.*", "", contig_id)
source_ind <- sub(":.*", "", genome_metadata$genome) %>% setNames(genome_metadata$genome)
idx_by_mag <- split(seq_along(mag_of), mag_of)
idx_by_mag <- idx_by_mag[names(idx_by_mag) %in% genome_metadata$genome]

min_reads <- 500; min_contigs <- 10; min_contig_reads <- 20
# Excess spread (log2 units) of the contig-level abundance ratio between two samples,
# after subtracting the spread expected from Poisson counting noise. ~0 if all contigs of
# the bin change by the same factor (single population); larger if contigs move
# independently (mixed populations). Independent of sequencing depth.
disp_index <- function(a, b) {
  k <- (a + b) >= min_contig_reads
  if (sum(k) < min_contigs) return(NA_real_)
  a <- a[k]; b <- b[k]
  lr <- log2((a + .5) / (b + .5))
  noise <- mean(1/(a + .5) + 1/(b + .5)) / log(2)^2
  sqrt(max(0, var(lr) - noise))
}
pair_disp <- function(rows, samples) {
  s <- samples[colSums(mat[rows, samples, drop = FALSE]) >= min_reads]
  if (length(s) < 2) return(NULL)
  cb <- combn(s, 2)
  tibble(s1 = cb[1,], s2 = cb[2,],
         disp = map2_dbl(cb[1,], cb[2,], ~ disp_index(mat[rows, .x], mat[rows, .y])))
}
ind_samples <- split(smeta$Tube_code, smeta$individual)

real <- map_dfr(names(idx_by_mag), function(g) {
  r <- pair_disp(idx_by_mag[[g]], intersect(ind_samples[[source_ind[g]]], colnames(mat)))
  if (is.null(r)) return(NULL)
  mutate(r, genome = g, individual = source_ind[g])
}) %>% filter(!is.na(disp))

# artificial chimeras: pool contigs of two different MAGs from the same individual
chim <- map_dfr(split(names(idx_by_mag), source_ind[names(idx_by_mag)]), function(gs) {
  if (length(gs) < 2) return(NULL)
  pairs <- combn(gs, 2, simplify = FALSE)
  pairs <- pairs[sample(length(pairs), min(length(pairs), 40))]
  map_dfr(pairs, function(pr) {
    rows <- c(idx_by_mag[[pr[1]]], idx_by_mag[[pr[2]]])
    r <- pair_disp(rows, intersect(ind_samples[[source_ind[pr[1]]]], colnames(mat)))
    if (is.null(r)) return(NULL)
    mutate(r, genome = paste(pr, collapse = "+"), individual = source_ind[pr[1]])
  })
}) %>% filter(!is.na(disp))

per_mag <- real %>% group_by(genome, individual) %>%
  summarise(n_pairs = n(), median_disp = median(disp), max_disp = max(disp), .groups = "drop")
per_chim <- chim %>% group_by(genome, individual) %>%
  summarise(n_pairs = n(), median_disp = median(disp), max_disp = max(disp), .groups = "drop")
thr <- quantile(per_chim$median_disp, 0.05)       # 95% of artificial chimeras fall above
per_mag <- per_mag %>% mutate(below_chimera_5pct = median_disp < thr)

# ---- 2. host association of MAGs ---------------------------------------------
gc_mat <- genome_counts %>% column_to_rownames("genome") %>% as.matrix()
gc_mat <- gc_mat[, colnames(gc_mat) %in% smeta$Tube_code]
rel <- sweep(gc_mat, 2, colSums(gc_mat), "/")
smeta_g <- smeta %>% filter(Tube_code %in% colnames(rel))
assoc <- map_dfr(rownames(rel), function(g) {
  if (!any(smeta_g$individual == source_ind[g])) return(NULL)   # source individual without samples
  own <- smeta_g$Tube_code[smeta_g$individual == source_ind[g]]
  oth <- smeta_g$Tube_code[smeta_g$individual != source_ind[g]]
  tibble(genome = g, individual = source_ind[g],
         own_mean_rel = mean(rel[g, own]), other_mean_rel = mean(rel[g, oth]),
         own_prevalence = mean(rel[g, own] > 0), other_prevalence = mean(rel[g, oth] > 0),
         top_sample_is_own = smeta_g$individual[which.max(rel[g, smeta_g$Tube_code])] == source_ind[g],
         n_timepoints_detected_own = sum(rel[g, own] > 0))
}) %>% mutate(own_to_other_ratio = own_mean_rel / other_mean_rel)

# ---- 3. samples cluster by individual across time points ---------------------
bc <- as.matrix(vegdist(t(rel), method = "bray"))
sm <- smeta_g[match(colnames(bc), smeta_g$Tube_code), ]
pairs_d <- which(upper.tri(bc), arr.ind = TRUE) %>% as_tibble() %>%
  mutate(d = bc[cbind(row, col)],
         same_ind = sm$individual[row] == sm$individual[col],
         same_tp = sm$time_point[row] == sm$time_point[col],
         same_species = sm$species[row] == sm$species[col])
dist_summary <- pairs_d %>% filter(same_species) %>%
  group_by(same_ind, same_tp) %>% summarise(n = n(), mean_bc = mean(d), median_bc = median(d), .groups = "drop")
nn <- map_dfr(seq_len(ncol(bc)), function(i) {
  o <- setdiff(seq_len(ncol(bc)), i); o <- o[sm$time_point[o] != sm$time_point[i]]  # other time points only
  j <- o[which.min(bc[i, o])]
  tibble(sample = colnames(bc)[i], time_point = sm$time_point[i], nn_same_individual = sm$individual[j] == sm$individual[i])
})
nn_chance <- mean(sapply(seq_len(ncol(bc)), function(i) {
  o <- setdiff(seq_len(ncol(bc)), i); o <- o[sm$time_point[o] != sm$time_point[i]]
  mean(sm$individual[o] == sm$individual[i])}))
set.seed(123)
perm <- adonis2(as.dist(bc) ~ individual, data = sm, permutations = 9999)
set.seed(123)
perm_after <- { k <- sm$time_point %in% c("FMT1","FMT2"); adonis2(as.dist(bc[k,k]) ~ individual, data = sm[k,], permutations = 9999) }

# matched-time-point test: for each ordered pair of time points, is a sample closer to its own
# individual's sample than to other same-species individuals' samples? (percentile rank of own)
rank_tbl <- map_dfr(seq_len(ncol(bc)), function(i) {
  map_dfr(setdiff(levels(sm$time_point), as.character(sm$time_point[i])), function(tp) {
    o <- which(sm$time_point == tp & sm$species == sm$species[i])
    own <- o[sm$individual[o] == sm$individual[i]]
    if (length(own) != 1 || length(o) < 3) return(NULL)
    tibble(sample = colnames(bc)[i], from = as.character(sm$time_point[i]), to = tp,
           pct_closer = mean(bc[i, o[o != own]] > bc[i, own]))   # share of other individuals farther than own
  })
})
rank_summary <- rank_tbl %>% summarise(n = n(), mean_pct = mean(pct_closer), median_pct = median(pct_closer), share_closest = mean(pct_closer == 1))
rank_wilcox <- wilcox.test(rank_tbl$pct_closer - 0.5)$p.value
# the 824 comparisons are not independent (several per sample and per individual), so the
# reported test averages them per individual (n = 27) and uses a signed-rank test against 0.5
rank_ind <- rank_tbl %>% left_join(smeta %>% select(Tube_code, individual), by = c("sample" = "Tube_code")) %>%
  group_by(individual) %>% summarise(mean_pct = mean(pct_closer), .groups = "drop")
rank_ind_wilcox <- wilcox.test(rank_ind$mean_pct - 0.5)$p.value
write_tsv(rank_ind, file.path(out, "own_vs_others_rank_per_individual.tsv"))
write_tsv(rank_tbl, file.path(out, "own_vs_others_rank.tsv"))

# ---- outputs -----------------------------------------------------------------
write_tsv(per_mag, file.path(out, "contig_dispersion_per_mag.tsv"))
write_tsv(per_chim, file.path(out, "contig_dispersion_artificial_chimeras.tsv"))
write_tsv(assoc, file.path(out, "host_association_per_mag.tsv"))
write_tsv(dist_summary, file.path(out, "bray_curtis_within_between.tsv"))
write_tsv(nn, file.path(out, "nearest_neighbour_per_sample.tsv"))

summ <- c(
  sprintf("MAGs in catalogue: %d; MAGs with >=2 eligible samples in source individual: %d", nrow(genome_metadata), nrow(per_mag)),
  sprintf("Sample pairs evaluated: real %d; artificial chimeras %d (%d chimeric bins)", nrow(real), nrow(chim), nrow(per_chim)),
  sprintf("Dispersion index, real MAGs: median %.2f (IQR %.2f-%.2f)", median(per_mag$median_disp), quantile(per_mag$median_disp,.25), quantile(per_mag$median_disp,.75)),
  sprintf("Dispersion index, artificial chimeras: median %.2f (IQR %.2f-%.2f); 5th percentile %.2f", median(per_chim$median_disp), quantile(per_chim$median_disp,.25), quantile(per_chim$median_disp,.75), thr),
  sprintf("Real MAGs below chimera 5th percentile: %d / %d (%.1f%%)", sum(per_mag$below_chimera_5pct), nrow(per_mag), 100*mean(per_mag$below_chimera_5pct)),
  sprintf("Real MAGs: share of max-pair dispersion above chimera median: %d", sum(per_mag$max_disp > median(per_chim$median_disp))),
  sprintf("Host association: median own/other mean rel. abundance ratio %.1f; MAGs with ratio>1: %d / %d", median(assoc$own_to_other_ratio[is.finite(assoc$own_to_other_ratio)]), sum(assoc$own_to_other_ratio > 1, na.rm=TRUE), nrow(assoc)),
  sprintf("Host association: highest-abundance sample belongs to source individual for %d / %d MAGs (%.1f%%)", sum(assoc$top_sample_is_own), nrow(assoc), 100*mean(assoc$top_sample_is_own)),
  sprintf("Prevalence in own individual's samples median %.2f vs other individuals %.2f", median(assoc$own_prevalence, na.rm=TRUE), median(assoc$other_prevalence)),
  sprintf("Nearest-neighbour (other time point) is same individual: %.1f%% (chance %.1f%%)", 100*mean(nn$nn_same_individual), 100*nn_chance),
  sprintf("Matched time points: own-individual sample is closer than %.1f%% of other individuals on average (chance 50%%; n = %d comparisons, not independent: see the per-individual test below; comparison-level Wilcoxon p = %.2g); closest of all in %.1f%%", 100*rank_summary$mean_pct, rank_summary$n, rank_wilcox, 100*rank_summary$share_closest),
  sprintf("Matched time points, per individual (n = %d): mean %.1f%%; %d individuals above 50%%; Wilcoxon signed-rank p = %.2g", nrow(rank_ind), 100*mean(rank_ind$mean_pct), sum(rank_ind$mean_pct > 0.5), rank_ind_wilcox),
  sprintf("PERMANOVA individual (all samples): R2 = %.3f, p = %.3f; (R2 inflated by 27 levels) FMT1+FMT2 only: R2 = %.3f, p = %.3f", perm$R2[1], perm$`Pr(>F)`[1], perm_after$R2[1], perm_after$`Pr(>F)`[1]))
writeLines(summ, file.path(out, "summary.txt"))
cat(paste(summ, collapse = "\n"), "\n")
print(dist_summary)
print(table(nn$time_point, nn$nn_same_individual))

# ---- figure ------------------------------------------------------------------
p1 <- bind_rows(mutate(per_mag, set = "Real MAGs"), mutate(per_chim, set = "Artificial chimeras")) %>%
  ggplot(aes(median_disp, fill = set)) + geom_histogram(bins = 40, position = "identity", alpha = .6) +
  scale_fill_manual(values = c("Real MAGs" = "#008080", "Artificial chimeras" = "#d57d2c")) +
  labs(tag = "a", x = "Excess contig-ratio spread (log2)", y = "Bins", fill = NULL) + theme_minimal()
p2 <- assoc %>% filter(own_mean_rel > 0 | other_mean_rel > 0) %>%
  ggplot(aes(other_mean_rel + 1e-6, own_mean_rel + 1e-6)) + geom_point(alpha = .5, colour = "#008080") +
  geom_abline(linetype = 2) + scale_x_log10() + scale_y_log10() +
  labs(tag = "b", x = "Mean rel. abundance in other individuals", y = "Mean rel. abundance in source individual") + theme_minimal()
p3 <- pairs_d %>% filter(same_species) %>%
  mutate(grp = case_when(same_ind ~ "Same individual,\ndifferent time point", same_tp ~ "Different individual,\nsame time point", TRUE ~ "Different individual,\ndifferent time point")) %>%
  ggplot(aes(grp, d)) + geom_boxplot(fill = "grey90", outlier.size = .3) +
  labs(tag = "c", x = NULL, y = "Bray-Curtis dissimilarity") + theme_minimal()
ggsave(file.path(out, "mag_individual_association.pdf"), gridExtra::grid.arrange(p1, p2, p3, nrow = 1), width = 15, height = 4.5)
