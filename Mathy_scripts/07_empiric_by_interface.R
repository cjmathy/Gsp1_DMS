#!/usr/bin/env RScript

# run from command line using 'Rscript 07_empiric_by_interface.R' from
# the scripts directory of the project

## Author: Christopher Mathy
## Date: 2020-02-11
## Email: cjmathy@gmail.com
## Email: chris.mathy@ucsf.edu
## Description:
##  This script analyzes the empiric data by looking into partner interface information

##### PREPARE DATASETS #####
# load modules
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)
library(factoextra)

# set ordering of amino acids
amino_acid_order <- c('*','K','R','H','D','E','N',
                      'Q','S','T','P','G','C','A',
                      'V','I','L','M','F','Y','W')

# for these plots, we'll only work with the 15 generation selection
secstruct <-
  readr::read_delim('data/ss_df_wide.csv', delim=',', col_types=cols()) %>%
  select(aa_position, gsp1_aa, consensus_ss)

nucleotide_binding_site <- data.frame(
  'partner' = c(rep('P-Loop', 7), rep('Switch I', 7), rep('Switch II', 9)),
  'yeastresnum' = c(seq(17, 23), seq(39, 45), seq(67, 75)),
  'in_interface' = c(rep(T, 7), rep(T, 9), rep(T, 7)),
  stringsAsFactors=F)

interfaces <-
  read_delim('data/SASA_interfaces.txt', delim='\t', col_types=cols()) %>%
  mutate(deltarASA = ifelse(is.na(deltarASA), 0, deltarASA)) %>%
  mutate(in_interface = ifelse(deltarASA > 0.25, T, F)) %>%
  select(-deltarASA, -interface) %>%
  unique() %>%
  bind_rows(nucleotide_binding_site) %>%
  right_join(data.frame('yeastresnum' = seq(from = 1, to = 220)),
             by = 'yeastresnum') %>%
  complete(partner, yeastresnum) %>%
  filter(!is.na(partner)) %>%
  mutate(in_interface = ifelse(is.na(in_interface), F, T))



empiric <-
  readr::read_delim('data/empiric_data_cleaned.txt', delim='\t', col_types=cols()) %>%
  left_join(secstruct, by = c('position' = 'aa_position', 'aa_from' = 'gsp1_aa')) %>%
  left_join(interfaces, by = c('position' = 'yeastresnum')) %>%
  mutate(aa_from = factor(aa_from, levels = amino_acid_order),
         aa_to = factor(aa_to, levels = amino_acid_order),
         dataset = factor(dataset, levels = c('6gen', '15gen')))

gsp1_seq <- select(empiric, aa_from, position) %>% unique() %>% pull(aa_from) %>% as.character()

##### HEATMAP OF DATASET WITH INTERFACE INFORMATION #####

make_mat <- function(ds) {
  empiric %>%
    filter(dataset == ds) %>%
    select(position, aa_to, fitness_score) %>%
    unique() %>%
    pivot_wider(id_cols = aa_to, names_from = position, values_from = fitness_score) %>%
    arrange(desc(aa_to)) %>%
    column_to_rownames('aa_to') %>%
    as.matrix()
}

make_heatmap_WT_labeled <- function(mat, name, col) {
  Heatmap(mat=mat, name=name, col=col, show_heatmap_legend = F, na_col='black',
    # rows
    row_title = 'Gsp1 residue', cluster_rows = F, row_names_side = 'left',
    row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
    row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
    # column
    column_title = 'Gsp1 position', cluster_columns = F,
    column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
    column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
    column_names_side = 'top',
    column_labels = replace(
      colnames(mat), seq(from = 2, to = length(colnames(mat)), by = 2), ''),
    # # uncomment to include borders around WT
    # cell_fun = function(j, i, x, y, width, height, fill) {
    #   if ((rownames(mat)[i] == gsp1_seq[j])) {
    #       grid.rect(x, y, width, height, gp = gpar(col = 'gray', lwd = 2, fill = "transparent"))
    #   }
    # }
  )
}

col_fn <- colorRamp2(seq(-4, 4, by=2), RColorBrewer::brewer.pal(5, 'RdBu'))


mat_6gen <- make_mat('6gen')
mat_15gen <- make_mat('15gen')

hm1 <- make_heatmap_WT_labeled(
  mat = mat_6gen, name='6 gen selection',
  col = col_fn)
hm2 <- make_heatmap_WT_labeled(
  mat = mat_15gen, name='15 gen selection',
  col = col_fn)

# add mean score per position annotation
mean_annot_6gen = HeatmapAnnotation(
  `Mean 6gen` = colMeans(mat_6gen, na.rm = T, dims = 1), show_legend = F,
  col = list(`Mean 6gen` = col_fn),
  which = 'column', simple_anno_size = unit(2, 'mm'), annotation_name_side = 'left')
mean_annot_15gen = HeatmapAnnotation(
  `Mean 15gen` = colMeans(mat_15gen, na.rm = T, dims = 1), show_legend = F,
  col = list(`Mean 15gen` = col_fn),
  which = 'column', simple_anno_size = unit(2, 'mm'), annotation_name_side = 'left')

# make annotation of secondary structure
ss_annot = HeatmapAnnotation(
  `secondary structure` = empiric %>%
    select(position, consensus_ss) %>%
    unique() %>%
    column_to_rownames('position') %>%
    as.matrix(),
  show_annotation_name = F, simple_anno_size = unit(2, 'mm'),
  # NOTE for color: NA is also gray (should only be stop codon at position 220)
  col = list(`secondary structure` = c('L' = 'gray', 'H' = 'orange', 'E' = 'purple')),
  annotation_legend_param = list('nrow'=1)
)

# make annotation for interfaces

# first, we make an interface matrix
mat_iface <-
  empiric %>%
  select(position, partner, in_interface) %>%
  unique() %>%
  pivot_wider(names_from = partner, values_from = in_interface,
              values_fill = list(in_interface = FALSE)) %>%
  column_to_rownames('position') %>%
  as.matrix()

# next, we cluster the interface matrix and order it by its first PC
cormat_iface <- cor(mat_iface, use = "pairwise.complete.obs", method = "pearson")
dissim_iface <- as.dist((1 - cormat_iface)/2)
hc_iface <- hclust(dissim_iface, method = "average" )
pc1_iface <-
  cmdscale(dist(cormat_iface), eig = T, k = 1)$point %>%
  as_tibble(rownames = 'partner') %>%
  arrange(`V1`) %>% pull(partner)
hc_iface <- rotate(hc_iface, pc1_iface)
iface_order <- hc_iface$order

mat_iface <- mat_iface[,iface_order]

interface_annot = HeatmapAnnotation(
  `interface` = mat_iface,
  col = list(`interface` = c('TRUE' = 'black', 'FALSE' = 'white')),
  simple_anno_size = unit(4, 'mm'),
  annotation_name_side = 'left',
  border = T, show_legend = F,
  annotation_name_gp = gpar(fontsize = 8)
)

# Make the legend for all heatmaps
leg <- Legend(at = c(-4, -2, 0, 2, 4), direction = 'horizontal', title = 'EMPIRIC score',
              col_fun = col_fn)

# plot heatmaps with secondary structure a nnotation and interface annotation
pdf('plots/empiric_heatmaps_with_interfaces.pdf', width = 13, height = 10)
draw(hm1 %v% mean_annot_6gen %v% hm2 %v% mean_annot_15gen %v% ss_annot %v% interface_annot,
     heatmap_legend_list = leg, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

pdf('plots/empiric_15gen_with_interfaces.pdf', width = 13, height = 8)
draw(hm2 %v% mean_annot_15gen %v% ss_annot %v% interface_annot,
     heatmap_legend_list = leg, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()


##### HEATMAPS CUT BY INTERFACE #####

make_split_heatmaps <- function(partner, positions) {

  groups <- c('Non-interface positions', paste0(partner, ' interface'))
  split_vec <- rep(groups[1], 220)
  split_vec[positions] <- groups[2]

  mats <- list('6 generation\nselection' = mat_6gen,
               '15 generation\nselection' = mat_15gen)
  hms <- list()

  for (i in seq_along(mats)) {

    mat = mats[[i]]
    name = names(mats)[[i]]

    hms[[i]] <- Heatmap(mat=mat, name=name, show_heatmap_legend = F,
      col=col_fn, na_col = 'black',

      # rows
      row_title = 'Gsp1 residue', cluster_rows = F, row_names_side = 'left',
      row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
      row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),

      # columns
      cluster_columns = F, column_names_side = 'top',
      column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
      column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),

      # split
      column_split = factor(split_vec, levels = groups),
      cluster_column_slices = F, column_gap = unit(2, 'mm')
      )
  }
  draw(hms[[1]] %v% mean_annot_6gen %v% hms[[2]] %v% mean_annot_15gen %v% ss_annot %v% interface_annot,
       heatmap_legend_list = leg, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}

pdf('plots/interface_split_heatmaps.pdf', width = 18, height = 8)
interfaces %>%
  group_by(partner, in_interface) %>%
  summarise(positions = list(yeastresnum)) %>%
  filter(in_interface) %>%
  mutate(heatmaps = map2(partner, positions, make_split_heatmaps))
dev.off()


##### SINA PLOTS OF INTERFACES #####

plots <-
  empiric %>%
  select(mutant, fitness_score, partner, in_interface, dataset) %>%
  group_by(partner, dataset, in_interface) %>%
  mutate(n_NAs = sum(is.na(fitness_score)), n = n()) %>%
  ungroup() %>%
  filter(!is.na(fitness_score)) %>%
  group_by(partner) %>%
  nest() %>%
  mutate(plot = map2(
    data, partner, ~ ggplot(
      data = .x, aes(in_interface, fitness_score)) +
      ggforce::geom_sina() +
      facet_grid(cols = vars(dataset)) +
      ggtitle(paste0('Partner: ', partner)) +
      scale_x_discrete(breaks = c(T, F),
                       labels = c('Mutation in the interface',
                                  'Mutation NOT in the interface')) +
      xlab ('') + ylab('EMPIRIC Fitness Score')
    )
  )

# plot sina plots for both datasets, with one partner per page
pdf('plots/scores_by_interface_sina.pdf', width = 12, height = 8)
plots$plot
dev.off()

# include number of NAs as geom_text

##### CDF PLOTS OF INTERFACES #####
# TODO: KS?

plot_edcf <- function(dataset) {
  empiric %>%
    filter(dataset == dataset) %>%
    select(mutant, fitness_score, partner, in_interface) %>%
    group_by_at(vars(-fitness_score)) %>%
    filter(sum(is.na(fitness_score)) == 0) %>%
    ungroup() %>%
    ggplot(aes(x = fitness_score, color = in_interface)) +
    stat_ecdf(geom = "step") +
    facet_wrap(vars(partner)) +
    ggtitle(paste('CDFs of fitness scores by partner',
                  'for non-NA positions across both datasets',
                  paste0('dataset = ', dataset), sep = '\n')) +
    xlab('EMPIRIC Fitness Score') + ylab ('Cumulative Density') +
    theme_bw() + theme(legend.position="none")
}

p1 <- plot_edcf('6gen')
p2 <- plot_edcf('15gen')
p <- cowplot::plot_grid(p1, p2)

# include a plot to show which partners have the most NAs
p_partner_NA_barchart <-
  empiric %>%
  select(mutant, fitness_score, dataset, partner, in_interface) %>%
  group_by(partner, dataset, in_interface) %>%
  mutate(percent_NA = sum(is.na(fitness_score))/n()) %>%
  ungroup() %>%
  mutate(partner = factor(partner, levels = filter(., in_interface) %>%
                            arrange(desc(percent_NA)) %>%
                            pull(partner) %>% unique())) %>%
  select(partner, percent_NA, dataset, in_interface) %>%
  unique() %>%
  ggplot(aes(x = in_interface, y = percent_NA, fill = in_interface, group = dataset)) +
  geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
  scale_x_discrete(breaks = c(T, F), labels = c('T', 'F')) +
  xlab('Mutation in interface (split by dataset)') + ylab('percent of scores that are NA') +
  facet_grid(cols = vars(partner)) +
  theme_bw() + theme(legend.position = 'top')

pdf('plots/scores_by_interface_ecdf.pdf', width = 12, height = 12)
cowplot::plot_grid(p, p_partner_NA_barchart, nrow = 2, rel_heights = c(3,1))

# Remake the CDF plots, but only show the positive values (to look at beneficial mutations)
maxscore <- max(empiric$fitness_score, na.rm=T)
minscore <- min(empiric$fitness_score, na.rm=T)
cowplot::plot_grid(
  p1 + coord_cartesian(xlim=c(0, maxscore), ylim=c(0.4, 1)),
  p2 + coord_cartesian(xlim=c(0, maxscore), ylim=c(0.4, 1)) + theme(legend.position = 'right')
)

# You can also make a CDF just for negative values:
cowplot::plot_grid(p1 + xlim(0, maxscore), p2 + xlim(0, maxscore) + theme(legend.position = 'right'))
cowplot::plot_grid(p1 + xlim(1, maxscore), p2 + xlim(1, maxscore) + theme(legend.position = 'right'))

# Remake the CDF plots, but  only show the negative values (to look at deleterious mutations)
cowplot::plot_grid(
  p1 + coord_cartesian(xlim=c(minscore, 0), ylim=c(0, 0.5)),
  p2 + coord_cartesian(xlim=c(minscore, 0), ylim=c(0, 0.5)) + theme(legend.position = 'right')
)

# You can also make a CDF just of negative values:
cowplot::plot_grid(p1 + xlim(minscore, 0), p2 + xlim(minscore, 0) + theme(legend.position = 'right'))
cowplot::plot_grid(p1 + xlim(minscore, -1), p2 + xlim(minscore, -1) + theme(legend.position = 'right'))
dev.off()

