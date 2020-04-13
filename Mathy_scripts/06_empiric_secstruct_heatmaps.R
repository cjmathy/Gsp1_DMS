#!/usr/bin/env RScript

# run from command line using 'Rscript 06_empiric_secstruct_heatmaps.R' from
# the scripts directory of the project

## Author: Christopher Mathy
## Date: 2020-02-05
## Email: cjmathy@gmail.com
## Email: chris.mathy@ucsf.edu
## Description:
##  This script preprocesses structures of Ran GTPase downloaded
##  using the script 'download_data.sh'. It parses and aligns
##  PDBs using the packages bio3d

##### PREPARE DATASETS #####
# load modules
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

# set ordering of amino acids
amino_acid_order <- c('*','K','R','H','D','E','N',
                      'Q','S','T','P','G','C','A',
                      'V','I','L','M','F','Y','W')

# load in EMPIRIC data with the secondary structure annotations
secstruct <-
  readr::read_delim('data/ss_df_wide.csv', delim=',', col_types=cols()) %>%
  select(aa_position, gsp1_aa, consensus_ss)

empiric <-
  readr::read_delim('data/empiric_data_cleaned.txt', delim='\t', col_types=cols()) %>%
  left_join(secstruct, by = c('position' = 'aa_position', 'aa_from' = 'gsp1_aa')) %>%
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



##### HEATMAPS #####

# function for plotting a heatmap of empiric data
# input dataset should be a matrix, with 21 aa labels as rows
# and each residue position as columns
Heatmap_empiric <- function(mat, name, col, label_WT = F) {
  Heatmap(mat, name=name, col=col, show_heatmap_legend = F, na_col='black',
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
    # borders around WT
    cell_fun = function(j, i, x, y, width, height, fill) {

      if ((rownames(mat)[i] == gsp1_seq[j]) & label_WT) {
        grid.rect(x, y, width, height, gp = gpar(col = 'limegreen', lwd = 2, fill = "transparent"))
      }
    }
  )
}

# col_fn <- colorRamp2(c(-5,0,5), RColorBrewer::brewer.pal(3, 'BrBG'))
col_fn <- colorRamp2(seq(-4, 4, by=2), RColorBrewer::brewer.pal(5, 'RdBu'))


hms <- list()

# make heatmaps
hm_6gen <- Heatmap_empiric(mat_6gen, '6 generation\nselection (top)', col_fn)
hm_15gen <- Heatmap_empiric(mat_15gen, '15 generation\nselection (bottom)', col_fn)
hm_diff <- Heatmap_empiric(mat_diff, 'Diff. between\n6 and 15 gen', col_fn)
hm_6gen_lbld <- Heatmap_empiric(mat_6gen, '6 generation\nselection (top)', col_fn, label_WT=T)
hm_15gen_lbld <- Heatmap_empiric(mat_15gen, '15 generation\nselection (bottom)', col_fn, label_WT=T)
hm_diff_lbld <- Heatmap_empiric(mat_diff, 'Diff. between\n6 and 15 gen', col_fn, label_WT=T)

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

# make the EMPIRIC score legend
leg <- Legend(at = c(-4, -2, 0, 2, 4), direction = 'horizontal', title = 'EMPIRIC score', col_fun = col_fn)

# prepare heatmap lists
hm_lists <- list()
hm_lists[[1]] <- hm_6gen %v% mean_annot_6gen %v% hm_15gen %v% mean_annot_15gen %v% ss_annot
hm_lists[[2]] <- hm_6gen_lbld %v% mean_annot_6gen %v% hm_15gen_lbld %v% mean_annot_15gen %v% ss_annot
hm_lists[[3]] <- hm_6gen %v% mean_annot_6gen %v% hm_15gen %v% mean_annot_15gen %v% hm_diff %v% ss_annot
hm_lists[[4]] <- hm_6gen_lbld %v% mean_annot_6gen %v% hm_15gen_lbld %v% mean_annot_15gen %v% hm_diff_lbld %v% ss_annot

# plot heatmaps with secondary structure annotation
pdf('plots/empiric_heatmaps_with_secstruct.pdf', width = 13, height = 8)

for (i in seq_along(hm_lists)) {
  draw(hm_lists[[i]], heatmap_legend_list = leg, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
}
dev.off()


##### SCORE DISTRIBUTIONS WITH SECONDARY STRUCTURE INFO #####

# as a companion plot, show distributions of scores for each secondary structure
p1 <-
  empiric %>%
  filter(position != 220, dataset != 'diff') %>%
  select(fitness_score, consensus_ss, dataset) %>%
  ggplot(aes(x = fitness_score, fill = consensus_ss)) +
  geom_histogram(bins = 100, position='identity', color='black') +
  facet_grid(rows = vars(consensus_ss), cols = vars(dataset),
             labeller=labeller(consensus_ss = c('L' = 'Loop',
                                                'H' = 'Helix',
                                                'E' = 'Sheet'))) +
  xlab('Fitness score') + ylab('Count') +
  ggtitle('Distribution of fitness scores by secondary structure') +
  theme(strip.text.x = element_text(size=12), strip.text.y = element_text(size=12), legend.position="none") +
  geom_text(
    data = empiric %>%
      filter(dataset != 'diff', !is.na(consensus_ss)) %>%
      group_by(dataset, consensus_ss) %>%
      mutate(is_NA =  is.na(fitness_score),
             n_NA = sum(is.na(fitness_score)),
             n = n(),
             percent_NA = 100*sum(is.na(fitness_score))/n()) %>%
      select(consensus_ss, dataset, percent_NA, n) %>%
      unique(),
    aes(x = -8, y = 100, label = paste0(n, ' scores\n', round(percent_NA, 1), '% are NA')),
    vjust=-0.2, color="black"
  )

p2 <-
  empiric %>%
  filter(dataset != 'diff') %>%
  group_by(dataset) %>%
  mutate(is_NA =  is.na(fitness_score), n_NA = sum(is.na(fitness_score)), n = n()) %>%
  ungroup() %>%
  ggplot(aes(x = dataset, fill = is_NA)) +
  geom_bar(width = 0.5) +
  geom_text(data =
              empiric %>%
              filter(dataset != 'diff') %>%
              group_by(dataset) %>%
              mutate(is_NA =  is.na(fitness_score), n_NA = sum(is.na(fitness_score)), n = n()) %>%
              group_by(is_NA, add=T) %>%
              mutate(percent_NA = case_when(!is_NA ~ round((n-n_NA)/n, 3), is_NA ~ round(n_NA/n, 3)),
                     label_pos = case_when(is_NA ~ percent_NA*n, !is_NA ~ 1*n)) %>%
              select(dataset, is_NA, percent_NA, label_pos) %>%
              unique(),
            aes(y = label_pos, label = 100*percent_NA), vjust=-0.2, color="black") +
  ggtitle('# NAs in dataset') + xlab('Dataset') + ylab('Count') +
  scale_fill_manual(name = 'Fitness score', labels = c('Not NA', 'NA'), values = c('gray', 'black')) +
  theme(legend.position="bottom")

p <- cowplot::plot_grid(p1, p2, rel_widths = c(3,1))

ggsave(plot=p, height=6, width=10, filename = 'plots/score_hist_by_secstructure.pdf', useDingbats=FALSE)



##### SCORE SCATTERPLOTS BY SECONDARY STRUCTURE SPLIT #####

p <-
  empiric %>%
  filter(dataset != 'diff', !is.na(consensus_ss)) %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  ggplot(aes(x=`6gen`, y=`15gen`)) +
    geom_point() +
    facet_grid(cols = vars(consensus_ss),
               labeller=labeller(consensus_ss = c('L' = 'Loop', 'H' = 'Helix','E' = 'Sheet'))) +
    xlab('Fitness Scores, 6 generation selection') + ylab('Fitness Scores,\n15 generation selection')
ggsave(plot=p, height=3, width=9, filename = 'plots/scatterplots_by_secstructure.pdf', useDingbats=FALSE)


##### HEATMAP HIGHLIGHTING OUTLIERS #####

# make a mask by rowwise mutating "
outlier_mask <-
  empiric %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  mutate('is_outlier' = case_when((`6gen` < -7) & (!is.na(`15gen`)) ~ T, T ~ F)) %>%
  select(position, aa_to, is_outlier) %>%
  pivot_wider(id_cols = aa_to, names_from = position, values_from = is_outlier) %>%
  arrange(desc(aa_to)) %>%
  column_to_rownames('aa_to') %>%
  as.matrix()
# sum(outlier_mask)  # this should be 147, because there are 147 outliers

hm_diff_outliers <-
  Heatmap(mat_diff, name = 'score difference (6gen-15gen)', col = col_fn, na_col='black',

          # rows
          row_title = 'Gsp1 residue', cluster_rows = F, row_names_side = 'left',
          row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),

          # column
          column_title = 'Gsp1 position\nBoxes are outliers (6gen unfit, 15gen more fit)', cluster_columns = F,
          column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          column_names_side = 'top',
          column_labels = replace(colnames(mat_diff),
                                  seq(from = 2, to = length(colnames(mat_diff)), by = 2),
                                  ''),
          # legend
          heatmap_legend_param = list('legend_direction'='horizontal'),

          cell_fun = function(j, i, x, y, width, height, fill) {
            if (outlier_mask[i,j]) {
              grid.rect(x, y, width, height, gp = gpar(col = 'limegreen', lwd = 2, fill = "transparent"))
            }
          }
  )

pdf('plots/outlier_diff_heatmap.pdf', width = 13, height = 4.5)
draw(hm_diff_outliers %v% ss_annot, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
dev.off()
