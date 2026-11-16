#!/usr/bin/env Rscript
# Run from repository root. No network access; saved input objects only.
suppressPackageStartupMessages({library(tidyverse); library(ape)})
set.seed(20261115)
OUT <- 'manuscript_figures/figure2'
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
write <- function(x, name) write_tsv(as_tibble(x), file.path(OUT, paste0(name,'.tsv')))
a <- new.env(); load('data/data_03122025.Rdata', a)
b <- new.env(); load('data/objects_03122025.Rdata', b)
z <- new.env(); load('data/beta_16092025.Rdata', z)
m <- a$genome_metadata
stopifnot(!anyDuplicated(m$genome), !anyNA(m), all(m$phylum %in% names(a$phylum_colors)))
tree <- ladderize(keep.tip(a$genome_tree, intersect(a$genome_tree$tip.label,m$genome)))
# Equivalent to force.ultrametric(method='extend'): extend terminal branches only.
d <- node.depth.edgelength(tree); n <- Ntip(tree)
terminal <- which(tree$edge[,2] <= n)
tree$edge.length[terminal] <- tree$edge.length[terminal] + max(d[seq_len(n)]) - d[tree$edge[terminal,2]]
d <- node.depth.edgelength(tree)
ord <- tree$edge[tree$edge[,2]<=n,2]; pos <- integer(n); pos[ord] <- seq_len(n)-1
write(tibble(id=seq_len(n+tree$Nnode),depth=d,parent=tree$edge[match(seq_len(n+tree$Nnode),tree$edge[,2]),1],tip_pos=c(pos,rep(NA,tree$Nnode)),label=c(tree$tip.label,rep(NA,tree$Nnode))), 'tree_nodes')
write(m %>% mutate(novel=species=='s__'), 'mag_metadata')
write(tibble(phylum=names(a$phylum_colors),colour=unname(a$phylum_colors)), 'phylum_colors')
# Composition uses the requested data snapshot; species/type metadata use the Wild Rmd's final snapshot.
wild <- b$sample_metadata %>% filter(time_point=='Wild')
stopifnot(!anyDuplicated(wild$Tube_code),all(wild$Tube_code %in% names(a$genome_counts_filt)))
x <- as.matrix(a$genome_counts_filt[,wild$Tube_code]); ph <- m$phylum[match(a$genome_counts_filt$genome,m$genome)]
stopifnot(!anyNA(x),all(is.finite(x)),all(x>=0),all(colSums(x)>0),!anyNA(ph))
tss <- function(x) {
  stopifnot(!anyNA(x), all(is.finite(x)), all(x >= 0), all(colSums(x) > 0))
  sweep(x,2,colSums(x),'/')
}
# Toy check: two genomes, two samples -> columns (1/4,3/4), (1/2,1/2).
stopifnot(isTRUE(all.equal(tss(matrix(c(1,3,2,2),2)),matrix(c(.25,.75,.5,.5),2))))
rel <- tss(x); comp <- rowsum(rel,ph)

# Toy detection check: species groups of two samples and one sample.
toy <- matrix(c(1,0,0,1,1,0),nrow=2)
stopifnot(identical(as.numeric(rowMeans(toy>0)),c(2/3,1/3)),
          identical(as.numeric(rowMeans(toy[,1:2,drop=FALSE]>0)),c(.5,.5)))
# Species composition among the 27 wild samples uses
# detection frequencies (per-species denominators), normalized in the SVG.
prev <- tibble(genome=a$genome_counts_filt$genome)
for (sp in c('Podarcis_liolepis','Podarcis_muralis')) {
  ids_sp <- wild$Tube_code[wild$species==sp]
  stopifnot(length(ids_sp)>0)
  prev[[paste0('prev_',sp)]] <- 100*rowMeans(x[,ids_sp,drop=FALSE]>0)
}
write(m %>% mutate(novel=species=='s__') %>% left_join(prev,by='genome'), 'mag_metadata')
rownames(rel) <- a$genome_counts_filt$genome
write(as_tibble(rel,rownames='genome') %>%
  pivot_longer(-genome,names_to='sample',values_to='rel') %>% filter(rel>0) %>%
  left_join(m %>% select(genome,phylum),by='genome') %>% select(genome,sample,rel,phylum),
  'sample_mag_abundance')
write(m %>% count(phylum,name='n'), 'phylum_counts')
write(wild %>% count(species,name='n'), 'species_counts')

stopifnot(max(abs(colSums(comp)-1))<1e-12)
write(as_tibble(t(comp),rownames='sample') %>% pivot_longer(-sample,names_to='phylum',values_to='rel') %>% left_join(wild,by=c('sample'='Tube_code')), 'sample_composition')
# Neutral NMDS: same k, trymax and distances as beta_neutral_nmds; no ellipses in source.
distmat <- as.matrix(z$beta_q1n$S)
ids <- rownames(distmat)[rownames(distmat) %in% wild$Tube_code]
stopifnot(setequal(ids,wild$Tube_code))
fit <- vegan::metaMDS(as.dist(distmat[ids,ids]),trymax=500,k=2,verbosity=FALSE,trace=FALSE)
points <- as_tibble(vegan::scores(fit,display='sites'),rownames='sample') %>% left_join(wild,by=c('sample'='Tube_code')) %>% group_by(species) %>% mutate(x_cen=mean(NMDS1),y_cen=mean(NMDS2)) %>% ungroup()
stopifnot(!anyNA(points))
write(points,'nmds_points'); write(tibble(seed=20261115,stress=fit$stress,tries=fit$tries),'nmds_fit')
# Reproduce comunity_elem + commun_wilcox_elem_wild exactly; fail on missing values.
e <- b$element_gift %>% filter(time_point=='Wild') %>% select(-time_point,-type)
stopifnot(!anyNA(e))
e <- e %>% filter(rowSums(across(where(is.numeric),~ . !=0))>0) %>% select(Tube_code,where(~is.numeric(.) && sum(.)>0)) %>% left_join(b$sample_metadata[c(1,3)],by='Tube_code')
long <- e %>% pivot_longer(-c(Tube_code,species),names_to='trait',values_to='value')
stopifnot(!anyNA(long),n_distinct(long$species)==2)
tests <- long %>% group_by(trait) %>% summarise(p_value=wilcox.test(value~species)$p.value,.groups='drop') %>% mutate(p_adjust=p.adjust(p_value,method='BH'))
write(tests,'functional_tests')
means <- long %>% group_by(trait,species) %>% summarise(value=mean(value),.groups='drop') %>% pivot_wider(names_from=species,values_from=value)
fun <- means %>% inner_join(filter(tests,p_adjust<0.05),by='trait') %>% left_join(b$uniqueGIFT_db[c(1,3,4)],by=c('trait'='Code_element')) %>% arrange(Function) %>% mutate(Difference=Podarcis_muralis-Podarcis_liolepis,group_color=ifelse(Difference<0,'Podarcis_liolepis','Podarcis_muralis'))
stopifnot(nrow(fun)>0,!anyNA(fun))
write(fun,'functional_capacity'); write(long,'functional_values')
message(sprintf('EXPORT: %d tree tips; %d wild samples; composition sums %.12f..%.12f; NMDS stress %.8f; %d/%d BH-significant pathways',n,nrow(wild),min(colSums(comp)),max(colSums(comp)),fit$stress,nrow(fun),nrow(tests)))
