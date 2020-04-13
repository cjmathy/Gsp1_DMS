#!/usr/bin/env RScript

# run from command line using 'Rscript 08_empiric_vs_emap.R' from
# the scripts directory of the project

## Author: Christopher Mathy
## Date: 2020-02-14
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

interfaces <-
  read_delim('data/SASA_interfaces.txt', delim='\t', col_types=cols()) %>%
  select(-interface) %>%
  unique() %>%
  right_join(data.frame('yeastresnum' = seq(from = 1, to = 220)),
             by = 'yeastresnum') %>%
  complete(partner, yeastresnum) %>%
  filter(!is.na(partner)) %>%
  mutate(deltarASA = ifelse(is.na(deltarASA), 0, deltarASA)) # 0 means not in interface

emap <-
  read_tsv('data/gsp1_pEMAP_avg_merged_gene_names.txt', col_types = cols()) %>%
  gather(-Gene, key = strain, value = score) %>%
  separate(strain, 'library_gene', sep = ' - ', remove = TRUE, extra = 'drop') %>%
  separate(Gene, c('Gene', 'mutant'), sep = ' - ') %>%
  filter(! mutant %in% c('GSP1-NAT', 'T34N', 'NTER3XFLAG WT', 'CTER3XFLAG WT')) %>%
  group_by(mutant) %>%
  mutate(mean_emap_score = mean(score, na.rm = T),
         median_emap_score = median(score, na.rm = T),
         sd_emap_score = sd(score, na.rm = T),
         emap_score_range = max(score, na.rm = T) - min(score, na.rm = T),
         emap_lower_quantile = quantile(score, probs=0.1, na.rm = T),
         emap_upper_quantile = quantile(score, probs=0.9, na.rm = T)) %>%
  ungroup()

empiric <-
  readr::read_delim('data/empiric_data_cleaned.txt', delim='\t', col_types=cols()) %>%
  left_join(secstruct, by = c('position' = 'aa_position', 'aa_from' = 'gsp1_aa')) %>%
  left_join(interfaces, by = c('position' = 'yeastresnum')) %>%
  mutate(in_interface = deltarASA > 0) %>%
  mutate(aa_from = factor(aa_from, levels = amino_acid_order),
         aa_to = factor(aa_to, levels = amino_acid_order),
         dataset = factor(dataset, levels = c('6gen', '15gen'))) %>%
  mutate(in_interface = deltarASA > 0) %>%
  inner_join(select(emap, -Gene, -library_gene, -score) %>% unique(), by = 'mutant')

gsp1_seq <- select(secstruct, gsp1_aa, aa_position) %>% unique() %>% pull(gsp1_aa) %>% as.character()


##### SIMPLE EMAP METRICS #####

make_scatterplot <- function(var, legend.position) {
  ggplot(empiric, aes_string(x = var, y = 'fitness_score', color = 'dataset')) +
    geom_point() + theme(legend.position = legend.position)
}

vars <- c('mean_emap_score', 'median_emap_score','sd_emap_score',
          'emap_score_range', 'emap_lower_quantile', 'emap_upper_quantile')
plots <- list()

for (i in seq_along(vars)) {
  legend.position = ifelse(i == 6, 'bottom', 'none')
  plots[[i]] <- make_scatterplot(vars[[i]], legend.position)
}

p <- cowplot::plot_grid(plotlist = plots)
ggsave(plot=p, height=6, width=9, filename = 'plots/emap_metrics.pdf', useDingbats=FALSE)


##### WHICH GENES ARE CORRELATED WITH THE PROFILES #####

emap_mat <-
  emap %>%
  select(mutant, library_gene, score) %>%
  pivot_wider(names_from = library_gene, values_from = score) %>%
  column_to_rownames('mutant') %>%
  as.matrix()

empiric_profiles <-
  empiric %>%
  select(mutant, dataset, fitness_score) %>%
  unique() %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  arrange(mutant = factor(mutant, levels = rownames(emap_mat))) %>%
  column_to_rownames('mutant') %>%
  as.matrix()

# reasonable to compute a correlation because we have no NAs across 55 positions
nrow(empiric_profiles)
max(rowSums(is.na(empiric_profiles)))

# across 1444 genes, we have at least 93% non-NAs
max(rowSums(is.na(emap_mat))/ncol(emap_mat))

# so just remove columns which have NAs from the matrix, leaving 1221 genes
emap_mat_noNA <- emap_mat[, colSums(is.na(emap_mat)) == 0]

correlations <-
  cor(emap_mat_noNA, empiric_profiles, method = 'spearman') %>%
  as.data.frame() %>%
  rownames_to_column('gene')

p <-
  correlations %>%
  filter(abs(`6gen`) > 0.3 & abs(`15gen`) > 0.3) %>%
  pivot_longer(c(`6gen`, `15gen`), names_to = 'dataset', values_to = 'corr') %>%
  mutate(gene = fct_reorder(gene, corr, mean)) %>%
  ggplot(aes(x = gene, y = corr)) +
  geom_bar(stat = 'identity') +
  facet_grid(rows = vars(dataset)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(plot=p, height=6, width=20, filename = 'plots/correlation_genes.pdf', useDingbats=FALSE)


# now, what are their genetic interaction profiles like?
pos_genes <-
  correlations %>%
  filter(`6gen` > 0.3 & `15gen` > 0.3) %>%
  pivot_longer(c(`6gen`, `15gen`), names_to = 'dataset', values_to = 'corr') %>%
  mutate(gene = fct_reorder(gene, corr, mean)) %>%
  arrange(gene) %>%
  pull(gene) %>%
  unique()

neg_genes <-
  correlations %>%
  filter(`6gen` < -0.3 & `15gen` < -0.3) %>%
  pivot_longer(c(`6gen`, `15gen`), names_to = 'dataset', values_to = 'corr') %>%
  mutate(gene = fct_reorder(gene, corr, mean)) %>%
  arrange(gene) %>%
  pull(gene) %>%
  unique() %>%
  unlist()

pos_mat <- emap_mat[,colnames(emap_mat) %in% pos_genes]
neg_mat <- emap_mat[,colnames(emap_mat) %in% neg_genes]

# add empiric scores to this plot
empiric_mat <-
  empiric %>%
  select(mutant, dataset, fitness_score) %>%
  unique() %>%
  pivot_wider(names_from = 'dataset', values_from = 'fitness_score') %>%
  column_to_rownames('mutant') %>%
  as.matrix()

empiric_mat <- empiric_mat[row.names(pos_mat),]

# set colors for heatmap
GI_cyan <- '#0BC3E8'
GI_yellow <- '#FDFB00'
GI_black <- '#000000'

hm_empiric <-
  Heatmap(empiric_mat, name = '',
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          row_title = 'EMPIRIC scores', cluster_rows = T, row_names_side = 'left',
          row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_title = 'Yeast gene', cluster_columns = F, column_names_side = 'top',
          column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          show_heatmap_legend = F)
dev.off()



hm_pos <-
  Heatmap(pos_mat, name = 'EMAP of genes that correlate with EMPIRIC',
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          row_title = 'Gsp1 mutant', cluster_rows = T, row_names_side = 'left',
          row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_title = 'Yeast gene', cluster_columns = F, column_names_side = 'top',
          column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          heatmap_legend_param = list('legend_direction'='horizontal'))
hm_neg <-
  Heatmap(neg_mat,  name = 'EMAP of genes that anticorrelate with EMPIRIC',
          col = colorRamp2(c(-4, 0, 4), c(GI_cyan, GI_black, GI_yellow)),
          row_title = 'Gsp1 mutant', cluster_rows = T, row_names_side = 'left',
          row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_title = 'Yeast gene', cluster_columns = F, column_names_side = 'top',
          column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          heatmap_legend_param = list('legend_direction'='horizontal'))

pdf('plots/emap_of_similar_genes.pdf', height = 8, width = 8)
draw(hm_empiric + hm_pos, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
draw(hm_empiric + hm_neg, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()

