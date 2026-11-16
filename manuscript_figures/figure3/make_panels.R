# Figure 3 panels: (a) neutral alpha diversity, (b) dbRDA of neutral beta diversity.
# Run from the repository root (build_figure3.py does this):
#   Rscript manuscript_figures/figure3/make_panels.R
# Writes fig3a.svg / fig3b.svg (drawn at final size, mm) and audit TSVs into manuscript_figures/figure3/.
suppressMessages({library(tidyverse); library(vegan); library(nlme); library(emmeans); library(svglite)})
out <- "manuscript_figures/figure3"
load("data/data_03122025.Rdata"); load("data/alpha_16092025.Rdata"); load("data/beta_16092025.Rdata")

# colours as in Figure 1: warm-control orange, cold-intervention purple, cold-control teal
grp_lev <- c("Hot_control", "Treatment", "Control")
grp_lab <- c("Warm-control (WC)", "Cold-intervention (CI)", "Cold-control (CC)")
grp_col <- c(Hot_control = "#d57d2c", Treatment = "#7b3294", Control = "#008080")
tp_lev  <- c("Acclimation", "FMT1", "FMT2")
tp_shape <- c(Acclimation = 16, FMT1 = 17, FMT2 = 15)          # circle, triangle, square
tp_shape_fill <- c(Acclimation = 21, FMT1 = 24, FMT2 = 22)    # filled equivalents (centroids)
FONT <- "Helvetica"

meta <- sample_metadata %>% filter(time_point %in% tp_lev) %>%
  mutate(time_point = factor(time_point, levels = tp_lev), type = factor(type, levels = grp_lev))

## ---- panel a: neutral alpha diversity ------------------------------------------------------
am <- alpha_div %>% left_join(meta, by = join_by(sample == Tube_code)) %>% filter(!is.na(time_point))
stopifnot(nrow(am) == 79)
mod <- nlme::lme(neutral ~ type * time_point, random = ~ 1 | individual, data = am)   # as in 09_FMT.Rmd
ct <- as.data.frame(emmeans(mod, pairwise ~ time_point | type, adjust = "tukey")$contrasts) %>%
  mutate(star = case_when(p.value < .001 ~ "***", p.value < .01 ~ "**", p.value < .05 ~ "*", TRUE ~ "ns"))
write_tsv(ct, file.path(out, "alpha_contrasts.tsv"))
write_tsv(am %>% select(sample, individual, type, time_point, neutral), file.path(out, "alpha_values.tsv"))

lab_fac <- function(x) factor(grp_lab[match(as.character(x), grp_lev)], levels = grp_lab)
am <- am %>% mutate(group = lab_fac(type))
ymax <- max(am$neutral)
br <- ct %>% filter(star != "ns") %>%
  separate(contrast, c("t1", "t2"), sep = " - ") %>%
  mutate(x1 = match(t1, tp_lev), x2 = match(t2, tp_lev), span = abs(x2 - x1), group = lab_fac(type)) %>%
  group_by(type) %>% arrange(span, x1, .by_group = TRUE) %>%
  mutate(k = row_number(), y = ymax * (1.05 + 0.095 * (k - 1))) %>% ungroup()
ytop <- max(c(ymax * 1.12, br$y * 1.07))

pa <- ggplot(am, aes(time_point, neutral, colour = type, fill = type)) +
  geom_boxplot(alpha = .16, width = .58, outlier.shape = NA, linewidth = .35, fatten = 1.6) +
  geom_point(aes(shape = time_point), position = position_jitter(width = .13, height = 0, seed = 1), size = 1.15, alpha = .85, stroke = 0) +
  geom_segment(data = br, aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, linewidth = .3, colour = "#222222") +
  geom_segment(data = br, aes(x = x1, xend = x1, y = y - ymax * .02, yend = y), inherit.aes = FALSE, linewidth = .3, colour = "#222222") +
  geom_segment(data = br, aes(x = x2, xend = x2, y = y - ymax * .02, yend = y), inherit.aes = FALSE, linewidth = .3, colour = "#222222") +
  geom_text(data = br, aes(x = (x1 + x2) / 2, y = y + ymax * .005, label = star), inherit.aes = FALSE, size = 2.6, colour = "#222222", vjust = 0, family = FONT) +
  geom_text(data = tibble(group = factor(grp_lab, levels = grp_lab), type = factor(grp_lev, levels = grp_lev)),
            aes(x = 2, y = ytop * 1.065, label = group, colour = type), inherit.aes = FALSE, fontface = "bold", size = 2.75, family = FONT) +
  facet_wrap(~ group, nrow = 1) +
  scale_colour_manual(values = grp_col) + scale_fill_manual(values = grp_col) + scale_shape_manual(values = tp_shape) +
  scale_y_continuous(limits = c(0, ytop * 1.1), breaks = seq(0, 60, 10), expand = c(0, 0)) +
  labs(x = NULL, y = "Neutral alpha diversity (q = 1)") +
  theme_classic(base_size = 7, base_family = FONT) +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_blank(),
        axis.text = element_text(colour = "#222222", size = 6.5), axis.title = element_text(size = 7),
        panel.grid.major.y = element_line(colour = "#ebebeb", linewidth = .25), panel.spacing = unit(4, "mm"),
        axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3))
svglite(file.path(out, "fig3a.svg"), width = 180 / 25.4, height = 56 / 25.4)
print(pa); dev.off()

## ---- panel b: dbRDA of neutral beta diversity, centroid trajectories ---------------------------
keep <- meta$Tube_code
d <- as.matrix(beta_q1n$S); d <- as.dist(d[rownames(d) %in% keep, colnames(d) %in% keep])
mb <- meta %>% arrange(match(Tube_code, labels(d)))
cca <- capscale(d ~ time_point * type, data = mb)                         # as in 09_FMT.Rmd
imp <- summary(cca)$cont$importance
pc <- round(100 * imp["Proportion Explained", 1:2], 1)
sc <- as.data.frame(scores(cca, display = "sites")) %>% rownames_to_column("Tube_code") %>%
  left_join(mb, by = "Tube_code") %>% mutate(group = lab_fac(type))
cen <- sc %>% group_by(type, group, time_point) %>%
  summarise(CAP1 = mean(CAP1), CAP2 = mean(CAP2), n = n(), .groups = "drop")
write_tsv(sc %>% select(Tube_code, individual, type, time_point, CAP1, CAP2), file.path(out, "dbrda_points.tsv"))
write_tsv(cen, file.path(out, "dbrda_centroids.tsv"))
write_tsv(tibble(axis = c("CAP1", "CAP2"), pct_total_variance = as.numeric(pc)), file.path(out, "dbrda_axes.tsv"))

# arrows between consecutive centroids, shortened (in mm) so the heads sit next to the centroid markers
xr <- range(c(sc$CAP1)) + c(-.12, .12); yr <- range(sc$CAP2) + c(-.12, .12)
W_mm <- 136; H_mm <- 78                  # approximate drawing area of the plot (mm), used to shorten arrows
sx <- W_mm / diff(xr); sy <- H_mm / diff(yr)
# Curved connectors: a quadratic Bezier between the two centroid centres (bowing to the left of the direction of
# travel), computed in mm so that the trimming is circular. The curve is then cut where it leaves the start marker
# and where it reaches the end marker, so it starts at the edge of one marker and the arrow tip stops at the edge of
# the next, heading towards its centre.
R_START <- 2.7; R_END <- 3.1      # mm from the marker centre (marker radius ~2.1 mm + gap; the tip ends R_END away)
K <- 0.46                          # bow of the control point, as a fraction of the segment length
cc <- cen %>% arrange(type, time_point) %>% group_by(type, group) %>%
  mutate(x2 = lead(CAP1), y2 = lead(CAP2), to = lead(time_point)) %>% filter(!is.na(x2)) %>% ungroup()
conn <- pmap_dfr(list(seq_len(nrow(cc)), cc$CAP1, cc$CAP2, cc$x2, cc$y2, as.character(cc$type), as.character(cc$time_point)),
  function(i, x1, y1, x2, y2, ty, fr) {
    p0 <- c((x1 - xr[1]) * sx, (y1 - yr[1]) * sy); p2 <- c((x2 - xr[1]) * sx, (y2 - yr[1]) * sy)
    dv <- p2 - p0; len <- sqrt(sum(dv^2)); nrm <- c(-dv[2], dv[1]) / len
    p1 <- (p0 + p2) / 2 + nrm * K * len
    t <- seq(0, 1, length.out = 600)
    pts <- cbind((1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1],
                 (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2])
    d0 <- sqrt((pts[, 1] - p0[1])^2 + (pts[, 2] - p0[2])^2); d2 <- sqrt((pts[, 1] - p2[1])^2 + (pts[, 2] - p2[2])^2)
    i0 <- min(which(d0 >= R_START)); i1 <- max(which(d2 >= R_END))
    tibble(id = i, type = ty, from = fr, x = xr[1] + pts[i0:i1, 1] / sx, y = yr[1] + pts[i0:i1, 2] / sy)
  }) %>% mutate(type = factor(type, levels = grp_lev))
write_tsv(cc %>% select(type, from = time_point, to, CAP1, CAP2, x2, y2), file.path(out, "dbrda_trajectories.tsv"))

pb <- ggplot() +
  geom_hline(yintercept = 0, linewidth = .25, colour = "#c8c8c8", linetype = "dotted") +
  geom_vline(xintercept = 0, linewidth = .25, colour = "#c8c8c8", linetype = "dotted") +
  geom_point(data = sc, aes(CAP1, CAP2, fill = type, shape = time_point), colour = "white", size = 1.9, alpha = .5, stroke = .15) +
  geom_path(data = conn, aes(x, y, group = id, colour = type), linewidth = .85, lineend = "round", linejoin = "round",
            arrow = arrow(length = unit(2.1, "mm"), angle = 24, type = "closed"), show.legend = FALSE) +
  geom_point(data = cen, aes(CAP1, CAP2, fill = type, shape = time_point), colour = "white", size = 5, stroke = .6) +
  geom_text(data = cen %>% mutate(lab = match(as.character(time_point), tp_lev), CAP2 = CAP2 + ifelse(time_point == "FMT1", .105 / sy, 0)),
            aes(CAP1, CAP2, label = lab), colour = "white", fontface = "bold", size = 2.5, family = FONT) +
  annotate("text", x = xr[1] + .03, y = yr[1] + .09, hjust = 0, vjust = 0, size = 2.3, colour = "#666666", family = FONT,
           label = "Large symbols: group centroids (1 \u2192 2 \u2192 3); arrows: centroid trajectories") +
  scale_colour_manual(values = grp_col, guide = "none") +
  scale_fill_manual(values = grp_col, breaks = grp_lev, labels = grp_lab, name = NULL) +
  scale_shape_manual(values = tp_shape_fill, name = NULL, labels = c("1  Acclimation", "2  FMT1", "3  FMT2")) +
  coord_cartesian(xlim = xr, ylim = yr, expand = FALSE) +
  labs(x = sprintf("CAP1 (%.1f%% of total variation)", pc[1]), y = sprintf("CAP2 (%.1f%%)", pc[2])) +
  theme_classic(base_size = 7, base_family = FONT) +
  theme(legend.position = "right", legend.box = "vertical", legend.spacing.y = unit(2, "mm"),
        legend.key.height = unit(4, "mm"), legend.text = element_text(size = 7), legend.margin = margin(0, 0, 0, 0),
        axis.text = element_text(colour = "#222222", size = 6.5), axis.title = element_text(size = 7),
        axis.line = element_line(linewidth = .3), axis.ticks = element_line(linewidth = .3)) +
  guides(fill = guide_legend(order = 1, override.aes = list(shape = 21, colour = "white", size = 3.2, alpha = 1, stroke = .4)),
         shape = guide_legend(order = 2, override.aes = list(colour = "white", fill = "#555555", size = 3, alpha = 1, stroke = .4)))
svglite(file.path(out, "fig3b.svg"), width = 180 / 25.4, height = (H_mm + 12) / 25.4)
print(pb); dev.off()
cat(sprintf("PANELS OK: a n=%d, b n=%d; CAP1 %.1f%%, CAP2 %.1f%%\n", nrow(am), nrow(sc), pc[1], pc[2]))
