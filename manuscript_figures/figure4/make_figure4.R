# Figure 4: normalised oxygen consumption across time in each experimental group (style of Figure 3a).
# Run from the repository root (build_figure4.py does this):
#   Rscript manuscript_figures/figure4/make_figure4.R
suppressMessages({library(tidyverse); library(readxl); library(nlme); library(emmeans); library(svglite)})
out <- "manuscript_figures/figure4"

grp_lev <- c("WC", "CI", "CC")
grp_lab <- c("Warm-control (WC)", "Cold-intervention (CI)", "Cold-control (CC)")
grp_col <- c(WC = "#d57d2c", CI = "#7b3294", CC = "#008080")        # colour scheme of Figure 1
tp_lev  <- c("Acclimation", "FMT1", "FMT2")
tp_shape <- c(Acclimation = 16, FMT1 = 17, FMT2 = 15)
FONT <- "Helvetica"
SCALE <- 1e4                                                         # values shown x 10^-4 ml s^-1 g^-1

r <- read_excel("data/respirometry_resp2.xlsx") %>%
  mutate(individual = as.factor(individual), cage = as.factor(cage), time_point = as.factor(time_point),
         type = as.factor(type), QC_normalized = suppressWarnings(as.numeric(QC_normalized)), r2 = as.numeric(r2)) %>%
  filter(!is.na(QC_normalized))
# same exclusions as 10_Physiology.Rmd (low-quality measurements, r2 and weight outliers)
d <- r %>% filter(weight != 5.19, weight != 5.91, weight != 5.03, r2 != 0.808, QC_normalized != 0.000036, r2 != 0.728,
                  !cage %in% c("Cold_17", "Cold_3", "Cold_6", "Cold_9", "Hot_2", "Hot_8")) %>%
  mutate(type = factor(as.character(type), levels = grp_lev),
         time_point = factor(recode(as.character(time_point), "0" = "Acclimation", "1" = "FMT1", "2" = "FMT2"), levels = tp_lev),
         value = QC_normalized * SCALE)
stopifnot(nrow(d) == 56)

# within-group pairwise contrasts: LMM with individual as random effect (as in the Methods), Tukey-adjusted
ct <- map_dfr(grp_lev, function(g) {
  m <- nlme::lme(QC_normalized ~ time_point, random = ~ 1 | individual, data = filter(d, type == g))
  as.data.frame(emmeans(m, pairwise ~ time_point, adjust = "tukey")$contrasts) %>% mutate(type = g, .before = 1)
}) %>% mutate(star = case_when(p.value < .001 ~ "***", p.value < .01 ~ "**", p.value < .05 ~ "*", TRUE ~ "ns"))
write_tsv(ct, file.path(out, "oxygen_contrasts.tsv"))
write_tsv(d %>% select(individual, cage, type, time_point, weight, r2, QC_normalized), file.path(out, "oxygen_values.tsv"))
write_tsv(d %>% count(type, time_point), file.path(out, "oxygen_n.tsv"))

lab_fac <- function(x) factor(grp_lab[match(as.character(x), grp_lev)], levels = grp_lab)
d <- d %>% mutate(group = lab_fac(type))
ymax <- max(d$value)
br <- ct %>% filter(star != "ns") %>%
  separate(contrast, c("t1", "t2"), sep = " - ") %>%
  mutate(across(c(t1, t2), ~ sub("^time_point", "", .x)),
         x1 = match(t1, tp_lev), x2 = match(t2, tp_lev), span = abs(x2 - x1), group = lab_fac(type)) %>%
  group_by(type) %>% arrange(span, x1, .by_group = TRUE) %>%
  mutate(k = row_number(), y = ymax * (1.05 + 0.095 * (k - 1))) %>% ungroup()
ytop <- max(c(ymax * 1.12, br$y * 1.07))

p <- ggplot(d, aes(time_point, value, colour = type, fill = type)) +
  geom_boxplot(alpha = .16, width = .58, outlier.shape = NA, linewidth = .35) +
  geom_point(aes(shape = time_point), position = position_jitter(width = .13, height = 0, seed = 1), size = 1.15, alpha = .85, stroke = 0) +
  geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = .3, colour = "#222222") +
  geom_segment(data = br, aes(x = x1, xend = x1, y = y - ymax * .02, yend = y), inherit.aes = FALSE, linewidth = .3, colour = "#222222") +
  geom_segment(data = br, aes(x = x2, xend = x2, y = y - ymax * .02, yend = y), inherit.aes = FALSE, linewidth = .3, colour = "#222222") +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = y + ymax * .005, label = star), inherit.aes = FALSE, size = 2.6, colour = "#222222", vjust = 0, family = FONT) +
  geom_text(data = tibble(group = factor(grp_lab, levels = grp_lab), type = factor(grp_lev, levels = grp_lev)),
            aes(x = 2, y = ytop * 1.065, label = group, colour = type), inherit.aes = FALSE, fontface = "bold", size = 2.75, family = FONT) +
  facet_wrap(~ group, nrow = 1) +
  scale_colour_manual(values = grp_col) + scale_fill_manual(values = grp_col) + scale_shape_manual(values = tp_shape) +
  scale_y_continuous(limits = c(0, ytop * 1.1), breaks = 0:6, expand = c(0, 0)) +
  labs(x = NULL, y = expression("Normalised O"[2] * " consumption (" * 10^-4 * " ml s"^-1 * " g"^-1 * ")")) +
  theme_classic(base_size = 7, base_family = FONT) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_blank(),
        axis.text = element_text(colour = "#222222", size = 6.5), axis.title = element_text(size = 7),
        panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = .25), panel.spacing = unit(4, "mm"),
        axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3))
svglite(file.path(out, "figure4.svg"), width = 180 / 25.4, height = 56 / 25.4)
print(p); dev.off()
cat("FIGURE 4 OK; n =", nrow(d), "\n")
