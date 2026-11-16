# Supplementary Figures S1-S6. Run from the repository root:
#   Rscript manuscript_figures/supplementary/make_supplementary_figures.R
suppressPackageStartupMessages({
  library(tidyverse); library(ggpubr); library(ggsignif); library(readxl); library(janitor)
})
set.seed(2802)
out_dir <- "manuscript_figures/supplementary"
load("data/data_03122025.Rdata"); load("data/objects_03122025.Rdata")
load("data/alpha_16092025.Rdata"); load("data/beta_16092025.Rdata")

col_cold <- "#008080"; col_warm <- "#d57d2c"; col_ci <- "#7b3294"
fam <- "Helvetica"
th <- theme_classic(base_size = 8, base_family = fam) +
  theme(strip.background = element_blank(), plot.subtitle = element_text(face = "italic", size = 7.5),
        plot.title.position = "plot", legend.position = "none")
W <- 180 / 25.4  # inches
save_fig <- function(p, name, h_mm) {
  ggsave(file.path(out_dir, paste0(name, ".png")), p, width = 180, height = h_mm, units = "mm", dpi = 300, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), p, width = 180, height = h_mm, units = "mm",
         device = if (capabilities("cairo")) cairo_pdf else pdf)
}

# ---- generic directional bar plot --------------------------------------
# diff = mean(A) - mean(B); positive -> higher in A
bar_plot <- function(dt, A, B, colA, colB, lab_A, lab_B, subtitle = NULL, title = NULL) {
  dt <- dt %>% mutate(dir = ifelse(Difference > 0, "A", "B"),
                      lab = ifelse(is.na(Function), Elements, Function),
                      lab = paste0(lab, " (", Elements, ")"))
  lim <- max(abs(dt$Difference)) * 1.05
  n <- nrow(dt)
  ggplot(dt, aes(x = fct_reorder(lab, Difference), y = Difference, fill = dir)) +
    geom_col(width = 0.75) + geom_hline(yintercept = 0, linewidth = 0.3) +
    scale_fill_manual(values = c(A = colA, B = colB)) +
    scale_x_discrete(expand = expansion(add = c(0.6, 1.8))) +
    scale_y_continuous(limits = c(-lim, lim)) +
    annotate("text", x = n + 1.5, y = lim, label = paste("higher in", lab_A), hjust = 1, size = 2.9, fontface = "bold", colour = colA, family = fam) +
    annotate("text", x = n + 1.5, y = -lim, label = paste("higher in", lab_B), hjust = 0, size = 2.9, fontface = "bold", colour = colB, family = fam) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = paste0("Mean difference in trait value (", lab_A, " − ", lab_B, ")"), title = title, subtitle = subtitle) +
    th + theme(panel.grid.major.x = element_line(linewidth = 0.15, colour = "grey85"))
}
fig_h <- function(n) max(70, 30 + 5.2 * n)
mean_diff <- function(df, grp, A, B) {   # df: Tube_code + traits + grp column
  df %>% select(-Tube_code) %>% group_by(.data[[grp]]) %>% summarise(across(everything(), mean), .groups = "drop") %>%
    pivot_longer(-1, names_to = "Elements", values_to = "m") %>% pivot_wider(names_from = 1, values_from = m) %>%
    mutate(Difference = .data[[A]] - .data[[B]], mean_A = .data[[A]], mean_B = .data[[B]]) %>%
    left_join(uniqueGIFT_db[c(1, 3, 4)], by = join_by(Elements == Code_element)) %>% arrange(Function)
}
prep <- function(filt) element_gift %>% filter({{ filt }}) %>% select(-time_point, -type) %>%
  filter(rowSums(across(where(is.numeric), ~ . != 0)) > 0) %>% select(Tube_code, where(~ is.numeric(.) && sum(.) > 0))
wilcox_sig <- function(wide, grp_df, grp) {
  wide %>% left_join(grp_df, by = "Tube_code") %>%
    pivot_longer(-c(Tube_code, all_of(grp)), names_to = "trait", values_to = "value") %>%
    group_by(trait) %>% summarise(p_value = wilcox.test(as.formula(paste("value ~", grp)))$p.value, .groups = "drop") %>%
    mutate(p_adjust = p.adjust(p_value, "BH")) %>% filter(p_adjust < 0.05)
}
cat("sample_metadata type labels:", paste(unique(sample_metadata$type), collapse = ", "), "\n")

# ---- S1 alpha diversity at acclimation ------------------------------------
a <- alpha_div %>% pivot_longer(-sample, names_to = "metric", values_to = "value") %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code)) %>% filter(time_point == "Acclimation") %>%
  mutate(metric = factor(metric, levels = c("richness", "neutral", "phylogenetic"),
                         labels = c("Richness", "Neutral", "Phylogenetic")),
         species = factor(species, levels = c("Podarcis_muralis", "Podarcis_liolepis"),
                          labels = c("P. muralis", "P. liolepis")))
cat("S1 samples per species:\n"); print(table(a$species[a$metric == "Richness"]))
s1 <- ggplot(a, aes(species, value, colour = species, fill = species)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.7) +
  scale_colour_manual(values = c("P. muralis" = col_cold, "P. liolepis" = col_warm)) +
  scale_fill_manual(values = c("P. muralis" = col_cold, "P. liolepis" = col_warm)) +
  facet_wrap(~metric, scales = "free_y") +
  stat_compare_means(method = "t.test", comparisons = list(c("P. muralis", "P. liolepis")), label = "p.format", size = 2.8, tip.length = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL, y = "Alpha diversity", subtitle = "Welch two-sample t-test between species") +
  th + theme(axis.text.x = element_text(face = "italic"), legend.position = "none")
save_fig(s1, "figS1", 75)

# ---- S2 species at acclimation ---------------------------------------------
w <- prep(time_point == "Acclimation")
sig <- wilcox_sig(w, sample_metadata[c("Tube_code", "species")], "species")
cat("S2: acclimation samples:", nrow(w), "; per species:", paste(names(table(sample_metadata$species[match(w$Tube_code, sample_metadata$Tube_code)])), table(sample_metadata$species[match(w$Tube_code, sample_metadata$Tube_code)]), collapse = ", "), "\n")
dt2 <- w %>% select(Tube_code, all_of(sig$trait)) %>% left_join(sample_metadata[c("Tube_code", "species")], by = "Tube_code") %>%
  mean_diff("species", "Podarcis_muralis", "Podarcis_liolepis")
cat("S2 bars:", nrow(dt2), " higher in muralis:", sum(dt2$Difference > 0), " higher in liolepis:", sum(dt2$Difference < 0), "\n")
s2 <- bar_plot(dt2, "m", "l", col_cold, col_warm, "P. muralis (cold-adapted)", "P. liolepis (warm-adapted)",
               subtitle = "Wilcoxon rank-sum test per trait, Benjamini-Hochberg p < 0.05")
save_fig(s2, "figS2", fig_h(nrow(dt2)))

# ---- S3/S4 Acclimation vs FMT1 within group: PAIRED per-individual plots -----------
# Paired signed-rank test per trait (as in 09_FMT.Rmd chunk paired_trait_tests_acc_fmt1), BH across traits.
acc_fmt1 <- function(ty, col, tag, name, n_show = 12) {
  w <- prep(time_point %in% c("FMT1", "Acclimation") & type == ty)
  md <- sample_metadata[c("Tube_code", "time_point")]
  sig_un <- wilcox_sig(w, md, "time_point")$trait            # unpaired, BH < 0.05
  long <- w %>% left_join(sample_metadata[c("Tube_code", "time_point", "individual")], by = "Tube_code") %>%
    pivot_longer(-c(Tube_code, time_point, individual), names_to = "trait", values_to = "value")
  wd <- long %>% select(-Tube_code) %>% pivot_wider(names_from = time_point, values_from = value, values_fn = mean) %>%
    filter(!is.na(Acclimation), !is.na(FMT1))
  n_pairs <- n_distinct(wd$individual)
  pp <- wd %>% group_by(trait) %>%
    summarise(p = tryCatch(suppressWarnings(wilcox.test(Acclimation, FMT1, paired = TRUE, exact = FALSE)$p.value),
                           error = function(e) NA_real_), .groups = "drop") %>%
    mutate(p_adj = p.adjust(p, "BH"), unpaired_sig = trait %in% sig_un) %>%
    left_join(uniqueGIFT_db[c(1, 3, 4)] %>% distinct(Code_element, .keep_all = TRUE), by = join_by(trait == Code_element)) %>%
    arrange(desc(unpaired_sig), p)
  n_fill <- max(n_show, sum(pp$unpaired_sig))
  shown <- pp %>% arrange(desc(unpaired_sig), p) %>% slice_head(n = n_fill) %>% arrange(p, trait)
  shown <- shown %>% mutate(fn = ifelse(is.na(Function), trait, str_trunc(Function, 56)),
                            title = paste0(fn, " (", trait, ")", ifelse(unpaired_sig, " *", "")),
                            title = factor(title, levels = title),
                            ann = paste0("p = ", formatC(p, format = "f", digits = 3), "\nBH p = ", formatC(p_adj, format = "f", digits = 3)))
  cat("\n", name, ":", ty, " pairs =", n_pairs, " unpaired significant =", length(sig_un),
      " smallest paired BH p =", signif(min(pp$p_adj, na.rm = TRUE), 3), "\n")
  print(as.data.frame(shown[c("trait", "Function", "p", "p_adj", "unpaired_sig")]))
  d <- wd %>% inner_join(shown[c("trait", "title", "ann")], by = "trait") %>%
    pivot_longer(c(Acclimation, FMT1), names_to = "tp", values_to = "value") %>%
    mutate(tp = factor(tp, levels = c("Acclimation", "FMT1")))
  mm <- d %>% group_by(title, tp) %>% summarise(value = mean(value), .groups = "drop")
  an <- shown %>% select(title, ann)
  p <- ggplot(d, aes(tp, value)) +
    geom_line(aes(group = individual), colour = "grey65", linewidth = 0.3) +
    geom_point(colour = col, size = 1.2, alpha = 0.8) +
    geom_line(data = mm, aes(group = 1), colour = "black", linewidth = 0.8) +
    geom_point(data = mm, shape = 23, fill = col, colour = "black", size = 2.6, stroke = 0.5) +
    geom_text(data = an, aes(x = 1.5, y = Inf, label = ann), inherit.aes = FALSE, vjust = 1.15,
              size = 2.2, family = fam, lineheight = 0.9) +
    facet_wrap(~title, ncol = 3, scales = "free_y", labeller = label_wrap_gen(34)) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.35))) +
    labs(x = NULL, y = "Trait value", title = tag,
         subtitle = paste0("Paired samples of the same individuals (n = ", n_pairs, " pairs); signed-rank test;\nno trait significant after ",
                           "Benjamini-Hochberg correction (smallest adjusted p = ", formatC(min(pp$p_adj, na.rm = TRUE), format = "f", digits = 2), ")"),
         caption = "* significant in an unpaired rank-sum comparison that treats samples as independent") +
    th + theme(strip.text = element_text(size = 6.5, face = "bold"), plot.caption = element_text(hjust = 0, size = 7),
               plot.caption.position = "plot", panel.spacing = unit(3, "mm"))
  save_fig(p, name, 34 + 42 * ceiling(nrow(shown) / 3))
  p
}
s3 <- acc_fmt1("Cold_control", col_cold, "Cold-control", "figS3")
s4 <- acc_fmt1("Cold_intervention", col_ci, "Cold-intervention", "figS4")

# ---- S5 CC vs CI at FMT2 -----------------------------------------------------
w <- prep(type %in% c("Cold_control", "Cold_intervention") & time_point == "FMT2")
md <- sample_metadata[c("Tube_code", "type")]
sig <- wilcox_sig(w, md, "type")
ty <- md$type[match(w$Tube_code, md$Tube_code)]
cat("S5 samples Cold_control =", sum(ty == "Cold_control"), " Cold_intervention =", sum(ty == "Cold_intervention"), "\n")
dt5 <- w %>% select(Tube_code, all_of(sig$trait)) %>% left_join(md, by = "Tube_code") %>% mean_diff("type", "Cold_control", "Cold_intervention")
cat("S5 bars:", nrow(dt5), "\n"); print(as.data.frame(dt5[c("Elements", "Function", "mean_A", "mean_B", "Difference")]))
s5 <- bar_plot(dt5, "A", "B", col_cold, col_ci, "cold-control", "cold-intervention",
               subtitle = "FMT2; Wilcoxon rank-sum test per trait, Benjamini-Hochberg p < 0.05")
save_fig(s5, "figS5", fig_h(nrow(dt5)))

# ---- S6 oxygen consumption at acclimation ---------------------------------------
resp <- read_excel("data/respirometry_resp2.xlsx") %>%
  mutate(across(c(individual, cage, time_point, type, population), as.factor),
         QC_normalized = as.numeric(QC_normalized), r2 = as.numeric(r2)) %>% filter(!is.na(QC_normalized))
acc_excl <- c("AH1_2nd_8", "LI1_2nd_4", "LI1_2nd_3", "LI1_2nd_10")
acc <- resp %>% filter(time_point == "0", !individual %in% acc_excl) %>%
  mutate(population = factor(population, levels = c("Cold", "Warm"), labels = c("Cold-adapted", "Warm-adapted")))
nn <- table(acc$population); cat("S6 n per population:", paste(names(nn), nn, collapse = ", "), "\n")
wt <- wilcox.test(QC_normalized ~ population, data = acc); cat("S6 Wilcoxon W =", wt$statistic, " p =", format(wt$p.value, digits = 4), "\n")
s6 <- ggplot(acc, aes(population, QC_normalized, colour = population, fill = population)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.25, linewidth = 0.4) + geom_jitter(width = 0.12, size = 1.3, alpha = 0.8) +
  scale_colour_manual(values = c("Cold-adapted" = col_cold, "Warm-adapted" = col_warm)) +
  scale_fill_manual(values = c("Cold-adapted" = col_cold, "Warm-adapted" = col_warm)) +
  geom_signif(comparisons = list(c("Cold-adapted", "Warm-adapted")), test = "wilcox.test", map_signif_level = FALSE,
              textsize = 2.8, colour = "black", family = fam, tip_length = 0.02,
              annotations = paste0("Wilcoxon p = ", format(round(wt$p.value, 3), nsmall = 3))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = expression("Mass-specific O"[2]*" consumption (ml s"^-1*" g"^-1*")"),
       subtitle = paste0("Four individuals with low-quality measurements excluded (n = ", nn[1], " cold-adapted, ", nn[2], " warm-adapted)")) +
  th
save_fig(s6, "figS6", 85)
