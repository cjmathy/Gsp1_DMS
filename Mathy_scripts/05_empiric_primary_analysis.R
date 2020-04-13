#!/usr/bin/env RScript

# run from command line using 'Rscript 05_empiric_primary_analysis.R' from
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

empiric <-
  readr::read_delim('data/empiric_data_cleaned.txt', delim='\t', col_types=cols()) %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  mutate(diff = `15gen`-`6gen`) %>%
  pivot_longer(cols = c('6gen', '15gen', 'diff'), names_to = 'dataset', values_to = 'fitness_score') %>%
  mutate(aa_from = factor(aa_from, levels = amino_acid_order),
         aa_to = factor(aa_to, levels = amino_acid_order),
         dataset = factor(dataset, levels = c('6gen', '15gen', 'diff')))

gsp1_seq <- select(empiric, aa_from, position) %>% unique() %>% pull(aa_from) %>% as.character()


##### SCORE DISTRIBUTIONS COMPARING DATASETS #####

# A good first question to ask is: what is the expected noise due to the assay,
# which we measure here as spread of fitness scores where the sequence is wild
# type (i.e. the mutant plasmid encodes a WT sequence in the EMPIRIC selection)
# we include the mean and sd of that subset of points in these plots

p1 <-
  empiric %>%
  filter(dataset != 'diff') %>%
  mutate(is_WT_seq = ifelse(aa_to == aa_from, T, F)) %>%
  ggplot(aes(x = fitness_score, fill = is_WT_seq)) +
  geom_histogram(bins = 75, color = 'black', alpha = 0.2,
                 # aes(y = ..density..),  # uncomment to include density interpolation
                 position='identity') +
  # geom_density(color = 'black', alpha = 0.2, position='identity') +
  geom_text(data =
              empiric %>%
              filter(dataset != 'diff', aa_to == aa_from) %>%
              mutate(is_WT_seq = T) %>%  # already filtered for this in last line
              group_by(dataset) %>%
              mutate(mean = round(mean(fitness_score, na.rm = T), 2),
                     sd = round(sd(fitness_score, na.rm = T), 2)) %>%
              ungroup() %>%
              select(dataset, mean, sd, is_WT_seq) %>%
              unique(),
            aes(x = -12, y = 300, label = paste0('WT sequence points (blue):\n  Mean = ', mean, '\n  SD = ', sd)),
            hjust = 0) +
  xlab('Fitness score') + ylab('Count') + ggtitle('Distribution of fitness scores') +
  facet_grid(rows = vars(dataset)) +
  theme(strip.text.y = element_text(size=12), legend.position="none")

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

p <- cowplot::plot_grid(p1, p2, rel_widths = c(3,2))

ggsave(plot=p, height=6, width=8, filename = 'plots/score_hist.pdf', useDingbats=FALSE)


##### SCORE SCATTERPLOTS COMPARING DATASETS #####

# EMPIRIC scores, not normalized
p1 <-
  empiric %>%
  filter(dataset != 'diff') %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  filter(! (is.na(`6gen`) | is.na(`15gen`))) %>%
  ggplot(aes(x=`6gen`, y=`15gen`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color='blue') +
  geom_abline(intercept = 0, slope = 2.5, color='red') +
  ggtitle('blue line: 15gen = 6gen, red line: 15gen = 2.5*6gen') +
  xlab('Fitness Scores, 6 generation selection') + ylab('Fitness Scores,\n15 generation selection')

# EMPIRIC scores, normalized
p2 <-
  empiric %>%
  filter(dataset != 'diff') %>%
  mutate(fitness_score_norm = case_when(dataset == '6gen' ~ fitness_score/6,
                                        dataset == '15gen' ~ fitness_score/15)) %>%
  select(-fitness_score) %>%
  pivot_wider(names_from = dataset, values_from = fitness_score_norm) %>%
  filter(! (is.na(`6gen`) | is.na(`15gen`))) %>%
  ggplot(aes(x=`6gen`, y=`15gen`)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color='red') +
  ggtitle('red line: 15gen/15 = 6gen/6') +
  xlim(-2, 1) + ylim(-2, 1) +
  xlab('Fitness Scores norm. to n_gen, 6 generation selection') + ylab('Fitness Scores norm. to n_gen,\n15 generation selection')

p <- cowplot::plot_grid(p1, p2)

ggsave(plot=p, height=6, width=12, filename = 'plots/scatterplot_comparing_selections.pdf', useDingbats=FALSE)




##### EXAMINING OUTLIERS (UNFIT IN 6GEN BUT NOT 15GEN)

# What is the population of points seen in the scatterplot that are strongly
# unfit (score < -7) in the 6 generation selection, but are not as unfit in
# the 15 generation selection
# (This is unexpected because generally the trend we saw is for scores to gain
# in absolute value in the second, longer selection)

# The difference map is plotted with outlier positions boxed in the
# script 06_empiric_secstruct_heatmaps.R
# Also, the outliers can be seen as primarily in loops and helices using
# the scatterplots that are stratified by secondary structure

# doesn't seem like there's an particular preference for amino acid to / from
# as also seen in this clustered heatmap

outlier_mat <-
  empiric %>%
  pivot_wider(names_from = dataset, values_from = fitness_score) %>%
  filter(`6gen` < -7, !is.na(`15gen`)) %>%
  select(aa_from, aa_to) %>%
  group_by(aa_from, aa_to) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  unique() %>%
  pivot_wider(names_from = aa_to, values_from = n, values_fill = list(n = 0)) %>%
  right_join(data.frame('aa_from' = amino_acid_order)) %>%
  select(aa_from, amino_acid_order) %>%
  mutate_all(replace_na, 0) %>%  # Cys does not appear in subset, so it gets added back in
  column_to_rownames('aa_from') %>%
  as.matrix()

# unclustered
hm1 <- Heatmap(
  outlier_mat, name = 'Num. of outliers', col = RColorBrewer::brewer.pal(4, 'Blues'),
  row_title = 'Original Amino Acid', column_title = 'Mutated Amino Acid',
  row_names_side = 'left', column_names_side = 'top',
  cluster_rows = FALSE, cluster_columns = FALSE)

# clustered
hm2 <- Heatmap(
  outlier_mat, name = 'Num. of outliers\n(6gen score is\nmuch more deleterious\nthan 15gen score)',
  col = RColorBrewer::brewer.pal(4, 'Blues'),
  row_title = 'Original Amino Acid', column_title = 'Mutated Amino Acid',
  row_names_side = 'left', column_names_side = 'top',
  row_dend_width = unit(5, 'mm'), column_dend_height = unit(5, 'mm'))

pdf('plots/outlier_diff_heatmap.pdf', width = 7, height = 4.5)
draw(hm1)
draw(hm2)
dev.off()

