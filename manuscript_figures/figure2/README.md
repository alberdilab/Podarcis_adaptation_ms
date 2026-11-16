# Manuscript Figure 2

From the repository root:

```sh
Rscript manuscript_figures/figure2/export_data.R && python3 -I manuscript_figures/figure2/build_figure2.py
```

Requires the project's installed tidyverse, ape and vegan packages and Inkscape
at the application path in the builder. Python uses only its standard library.
The builder invokes Inkscape with GDK_BACKEND=x11 for headless macOS export.
The builder exports manuscript_figures/final/figure2.pdf (text outlined),
figure2.png (300 dpi), and figure2.tiff (LZW, RGB without alpha). It requires
ImageMagick `magick`; the bundled Inkscape Ghostscript renders a separate
figure2_pdf_preview.png for checking. All TSVs are reproducible outputs.
The pipeline and outputs live in manuscript_figures/figure2/, separate from
knitr outputs in Figures/.

## Data and statistical choices

- a/b: data/data_03122025.Rdata. Actual metadata fields are lowercase
  completeness and contamination. Novel species means species == 's__'.
  Tree is pruned to metadata and ladderized; terminal edges are extended to the
  greatest root-to-tip distance, as in the source force.ultrametric extend method.
  Display radii use a square-root depth scale, following the sled-dogs reference,
  to open the dense tips; this is a drawing transformation, not a new phylogeny.
  Branches have the phylum colour when all descendant tips share a phylum, grey
  otherwise. Largest contiguous phylum runs receive curved labels; small clades
  have leaders to labels placed without overlapping one another inside the fan.
  Pseudomonadota uses a complete straight label with a white halo and leader
  instead of the previously clipped curved text path. Curved label paths have
  sufficient length to display all their glyphs.
  Five rings, inside to outside: novel species (dark green / known light grey),
  genome quality, GC content, genome size, species composition.
  The separate phylum and completeness rings were removed. Completeness remains
  cheaply available as quality-bar height (50–100%); contamination is its red
  colour (0–10%). GC uses the saved `gc` fractions, unsmoothed line and area.
  Genome size bars use Mb, with unlabelled reference gridlines every 2 Mb
  (at 2, 4 and 6 Mb); the scale text is omitted for description in the caption.
  Species composition is each species' detection frequency (positive counts /
  number of wild samples of that species), normalized across the two species
  per MAG; this weights species equally despite n=9 versus n=18. A MAG absent
  from both species has a light-grey composition ring. The overall-prevalence
  ring, key, and unused prev_all export have been removed. Per-species detection
  frequencies remain necessary for the species-composition ring. The tree uses
  the freed radius; adjacent ring tracks have approximately 0.7–0.9 mm gaps.
- Wild filter: time_point == 'Wild'. Final species and experimental group labels
  come from objects_03122025.Rdata, matching load order in 06_Wild_populations.Rmd.
  b uses TSS of genome_counts_filt from the data snapshot and sums by phylum.
  Every wild sample must have positive total abundance. Zero individual MAGs
  remain zero; missing taxonomy, non-finite counts, missing samples or metadata
  cause export failure rather than silent omission. `sample_mag_abundance.tsv`
  has genome, sample, rel, phylum, with only nonzero MAG tiles; white 0.035 mm
  separators distinguish tiles. Within each species bars sort by descending
  Bacteroidota proportion (dominant across all wild samples), then sample code.
  The dominant phylum is stacked at the bottom; within phyla MAGs follow tree
  order. Each bar sums to one. Phylum labels drop only the p__ prefix.
- c: beta_q1n$S from beta_16092025.Rdata, subset in matrix order to wild IDs;
  vegan::metaMDS(k=2, trymax=500), seed 20261115. The original source does not
  save coordinates or set a seed, so exact historical orientation/coordinates
  cannot be recovered. This rerun preserves the distance input and model.
  Points and species centroid spokes match the source; no ellipses/hulls there.
  nmds_fit.tsv records stress, seed and tries. Legends and stress are added.
- d: reproduces the source's removal of zero-only rows and columns, per-element
  two-sided wilcox.test defaults, BH correction across all remaining pathways,
  adjusted p < 0.05 selection, and species arithmetic means. No transformation.
  Ties invoke the same normal approximation and continuity correction as the
  source; R warns for 18 tied tests. functional_tests.tsv includes all tests,
  functional_values.tsv all tested values, functional_capacity.tsv plotted means,
  differences and significance values. No stars are shown in the source.
  Function labels are ordered by median difference, matching fct_reorder defaults;
  labels now show the element (after the first underscore), followed by
  "biosynthesis" or "degradation" only when that word occurs in the function
  prefix. The statistical result remains 37/143 BH-significant pathways.
  A builder-only EXCLUDED_FUNCTIONS constant removes Cellular structure and
  Appendages from display: exactly Teichoic acid (S0103) and Pilus (S0202).
  Panel d shows 35 bars, with one-line labels and 2.4 mm row spacing (84 mm
  total). No rows are removed from functional_capacity.tsv or functional_tests.tsv.
  Bar fills match that prefix by exact Function name against data/gift_colors.tsv:
  12 displayed functions match, no missing names, spelling discrepancies or fallbacks.
  Unexpected names would print a warning and receive an unused fallback colour.
  functional_plot_labels.tsv is a presentation-only audit of the 35 shown bars,
  containing
  trait, original Function, function prefix, displayed label and colour; the
  existing statistical TSVs have no new columns or changed bytes.
- Species colours in b/c and d’s direction notes remain: orange P. liolepis
  (#d57d2c), teal P. muralis (#008080). Panel d shows muralis minus liolepis:
  positive means higher in muralis and negative means higher in liolepis.
  Bar fills now encode GIFT functions rather than species; species-coloured
  direction notes replace the species swatches in d. The x axis now extends
  left to -0.12; ticks use 0.04 increments and include zero. The chart expands
  leftwards instead of rescaling: zero stays at x=145.24268 mm and the linear
  scale remains 143.36616 mm/unit over both negative and positive values.
  Direction notes occupy two lines in the empty upper-left and lower-right
  corners, with white halos. The builder asserts that neither note box
  intersects any of the 35 bar rectangles. The plot top moves to y=122 mm.
- Panel b uses a 54.9 mm tall barplot and a compact three-column phylum legend
  directly below it. Single-line species headers include data-derived n and
  align at their ink top with the outer species-composition ring (12.1 mm).
  The plot spans x=109–175 mm, reducing the y-axis inset. The legend bottom and panel a caption align at a shared
  105 mm baseline. The rings legend block has been removed; short ring names
  remain at the open gap, and ring descriptions belong in the manuscript caption.
  The NMDS plot is now 40 mm tall with unchanged coordinates and centroid
  spokes. Stress = 0.134 is centred above the plot at x=45, y=124 mm.
  Species names are placed inside the plot beside their centroids, with a
  collision check against point symbols and every centroid spoke. The separate
  species key is removed; the compact shape key remains below the plot.
  NMDS2 is centred at y=147 mm on the 127–167 mm y axis and rotated -90°
  (bottom to top), at x=2.5 mm, left of the tick labels. NMDS1 is unchanged.
  The heading-free function legend is enclosed by a #f2f2f2 rounded rectangle
  (2 mm corner radius), spanning the full c column at x=2–86 mm (84 mm wide),
  just before the gutter to d at x=88 mm. It starts at y=186 mm below the shape
  legend. Twelve single-line entries are arranged 6+6 in the same order:
  first six biosynthesis functions left; Metallophore and Antibiotic biosynthesis
  followed by the four degradation functions right. The 12 displayed colours
  are a display-only override (DISPLAY_COLOURS in the builder; data/gift_colors.tsv
  is unchanged) that keeps the warm/cool scheme but separates similar hues: CIELAB
  dE of the closest pair rises from 10.9 to 27.8 (6.5 to 10.3 under simulated
  protan/deutan vision). Text is 6 pt with 2 mm inner padding and 0.6 mm swatch-to-text
  gaps. No labels wrap and no headings are shown. All eight biosynthesis colours
  are warm and four degradation colours cool; no HSV exceptions were found.
  Row 2 height is the greater of c (plot plus legends) and d (bars plus axes).

## Drawing changes and validation

178 mm wide (218.0 mm tall; reduced from 224.0 mm), rows a/b and c/d, Helvetica/Arial, all text at least 6 pt. Source
statistics and plotted values are unchanged, except for the explicit NMDS seed.
Labels, legend placement and panel aspect ratios adapt to the four-panel page.
Only the lizard species names P. liolepis and P. muralis are italicized using
inline SVG tspans (Helvetica/Arial italic): b headers, c cluster labels, d's
species direction lines, and both names in d's x-axis caption. Header counts,
arrows, other direction words, punctuation and surrounding axis text stay upright.
Phylum/taxonomic labels, MAG wording, GIFT functions and pathway names stay upright.
Collision-box estimates include 0.6 mm of italic overhang, and rendered bounds
are checked after export. The final PDF still outlines all text, including italics.
All 35 displayed pathway labels are present; the other two remain in the source TSVs. Export asserts unique IDs, no missing
plotted data, positive sample totals and composition sums within 1e-12 of one;
a two-sample toy matrix checks TSS against hand-calculated proportions, and
a three-sample toy matrix checks detection frequencies.

All five existing c/d TSVs are byte-identical to those before the redesign (SHA-256 checked).
Both the SVG preview and Ghostscript render of the outlined PDF were inspected.
Known warnings: the unchanged Wilcoxon tests warn about ties; headless Inkscape
may warn about GTK objects and its inaccessible extension log; exports succeed.

Retained one-time intermediate: redesign_top.py (source rewrite helper, not part
of the pipeline; do not rerun). TSVs, including functional_values.tsv, are audit
exports. figure2_pdf_preview.png is the independent PDF inspection render.

## Display label audit (35 shown bars)

The two display-excluded significant elements are S0103 (Teichoic acid) and S0202 (Pilus).

| Trait | Original Function | Displayed label | Colour |
|---|---|---|---|
| B0219 | Amino acid biosynthesis_GABA | GABA biosynthesis | #f79039 |
| B0213 | Amino acid biosynthesis_Proline | Proline biosynthesis | #f79039 |
| B0216 | Amino acid biosynthesis_Tryptophan | Tryptophan biosynthesis | #f79039 |
| D0503 | Amino acid degradation_Cysteine | Cysteine degradation | #b9dfc6 |
| D0518 | Amino acid degradation_GABA | GABA degradation | #b9dfc6 |
| D0508 | Amino acid degradation_Lysine | Lysine degradation | #b9dfc6 |
| D0517 | Amino acid degradation_Ornithine | Ornithine degradation | #b9dfc6 |
| B0302 | Amino acid derivative biosynthesis_Betaine | Betaine biosynthesis | #df5f26 |
| B0307 | Amino acid derivative biosynthesis_Spermidine | Spermidine biosynthesis | #df5f26 |
| B1012 | Antibiotic biosynthesis_Fosfomycin | Fosfomycin biosynthesis | #f7946f |
| B0803 | Aromatic compound biosynthesis_Chorismate | Chorismate biosynthesis | #d53f8c |
| B0801 | Aromatic compound biosynthesis_Salicylate | Salicylate biosynthesis | #d53f8c |
| B0902 | Metallophore biosynthesis_Aerobactin | Aerobactin biosynthesis | #ef4c78 |
| D0604 | Nitrogen compound degradation_GlcNAc | GlcNAc degradation | #96d2b1 |
| D0611 | Nitrogen compound degradation_Phenylethylamine | Phenylethylamine degradation | #96d2b1 |
| D0603 | Nitrogen compound degradation_Urate | Urate degradation | #96d2b1 |
| B0103 | Nucleic acid biosynthesis_UDP/UTP | UDP/UTP biosynthesis | #f9c145 |
| B0102 | Nucleic acid biosynthesis_Uridylic acid (UMP) | Uridylic acid (UMP) biosynthesis | #f9c145 |
| B0604 | Organic anion biosynthesis_L-lactate | L-lactate biosynthesis | #d51c58 |
| D0206 | Polysaccharide degradation_Alpha galactan | Alpha galactan degradation | #4a98c7 |
| D0212 | Polysaccharide degradation_Arabinan | Arabinan degradation | #4a98c7 |
| D0207 | Polysaccharide degradation_Beta-galactan | Beta-galactan degradation | #4a98c7 |
| D0204 | Polysaccharide degradation_Chitin | Chitin degradation | #4a98c7 |
| D0208 | Polysaccharide degradation_Mixed-Linkage glucans | Mixed-Linkage glucans degradation | #4a98c7 |
| D0213 | Polysaccharide degradation_Mucin | Mucin degradation | #4a98c7 |
| D0205 | Polysaccharide degradation_Pectin | Pectin degradation | #4a98c7 |
| D0203 | Polysaccharide degradation_Starch | Starch degradation | #4a98c7 |
| D0202 | Polysaccharide degradation_Xyloglucan | Xyloglucan degradation | #4a98c7 |
| D0304 | Sugar degradation_D-Arabinose | D-Arabinose degradation | #53b5db |
| D0309 | Sugar degradation_Galactose | Galactose degradation | #53b5db |
| D0307 | Sugar degradation_L-Fucose | L-Fucose degradation | #53b5db |
| D0310 | Sugar degradation_NeuAc | NeuAc degradation | #53b5db |
| D0302 | Sugar degradation_Sucrose | Sucrose degradation | #53b5db |
| B0708 | Vitamin biosynthesis_Cobalamin (B12) | Cobalamin (B12) biosynthesis | #a61f5e |
| B0705 | Vitamin biosynthesis_Pyridoxal-P (B6) | Pyridoxal-P (B6) biosynthesis | #a61f5e |
