#!/usr/bin/env RScript

# run from command line using 'Rscript 09_fitness_per_position.R' from
# the scripts directory of the project

## Author: Christopher Mathy
## Date: 2020-02-18
## Email: cjmathy@gmail.com
## Email: chris.mathy@ucsf.edu
## Description:
##  This script

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

# load in EMPIRIC data with the secondary structure annotations
# secstruct <-
#   readr::read_delim('data/ss_df_wide.csv', delim=',', col_types=cols()) %>%
#   select(aa_position, gsp1_aa, consensus_ss)

empiric <-
  readr::read_delim('data/empiric_data_cleaned.txt', delim='\t', col_types=cols()) %>%
  # left_join(secstruct, by = c('position' = 'aa_position', 'aa_from' = 'gsp1_aa')) %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  mutate(diff = `15gen`-`6gen`) %>%
  pivot_longer(cols = c('6gen', '15gen', 'diff'), names_to = 'dataset', values_to = 'fitness_score') %>%
  mutate(aa_from = factor(aa_from, levels = amino_acid_order),
         aa_to = factor(aa_to, levels = amino_acid_order),
         dataset = factor(dataset, levels = c('6gen', '15gen', 'diff')))

# prepare a sequence vector for WT Gsp1 so we can box cells in heatmaps that are WT in sequence (i.e. T34T)
gsp1_seq <- select(empiric, aa_from, position) %>% unique() %>% pull(aa_from) %>% as.character()

# convert dataframes to matrices for each dataset, for plotting
empiric_to_mat <- function(df, ds) {
  df %>%
    filter(dataset == ds) %>%
    select(position, aa_to, fitness_score) %>%
    pivot_wider(id_cols = aa_to, names_from = position, values_from = fitness_score) %>%
    arrange(desc(aa_to)) %>%
    column_to_rownames('aa_to') %>%
    as.matrix()
}
mat_6gen <- empiric_to_mat(empiric, '6gen')
mat_15gen <- empiric_to_mat(empiric, '15gen')
mat_diff <- empiric_to_mat(empiric, 'diff')







# cluster the 6gen heatmap
dissim <- get_dist(t(mat_6gen), method = "euclidean")
hc <- hclust(dissim, method = "ward.D2" )

# can rotate using first principle coordinate
pc1 <-
  cmdscale(dissim, eig = T, k = 1)$point %>%
  as_tibble(rownames = 'position') %>%
  arrange(`V1`) %>% pull(position)
hc <- rotate(hc, pc1)

# could also ratate using the mean score at each position
order_by_mean <- colMeans(mat_6gen, na.rm = T) %>% sort() %>% names()

order <- hc$order

# plot the heatmap
col_fn <- colorRamp2(seq(-4, 4, by=2), RColorBrewer::brewer.pal(5, 'RdBu'))
hm <-
  Heatmap(mat_15gen, name='clustered 15gen', col=col_fn, show_heatmap_legend = T, na_col = 'black',
        row_title = 'Gsp1 residue', cluster_rows = F, row_names_side = 'left',
        row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
        row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
        column_title = 'Gsp1 position',
        # cluster_columns = as.dendrogram(hc),
        column_order = order_by_mean,
        column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
        column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
        column_names_side = 'top',
        heatmap_legend_param = list('legend_direction' = 'horizontal')
        )

pdf('plots/ordered_heatmap.pdf', width = 16, height = 6)
draw(hm, heatmap_legend_side = "bottom")
dev.off()
