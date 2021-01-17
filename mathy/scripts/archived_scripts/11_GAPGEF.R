#!/usr/bin/env RScript

# run from command line using 'Rscript 11_GAPGEF.R' from
# the scripts directory of the project

## Author: Christopher Mathy
## Date: 2020-02-21
## Email: cjmathy@gmail.com
## Email: chris.mathy@ucsf.edu
## Description:
##  This script analyzes the empiric data by looking into GAP GEF interface  and kinetics

##### PREPARE DATASETS #####
# load modules
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)
library(factoextra)
library(scales)

# set ordering of amino acids
amino_acid_order <- c('*','K','R','H','D','E','N',
                      'Q','S','T','P','G','C','A',
                      'V','I','L','M','F','Y','W')

# for these plots, we'll only work with the 15 generation selection
secstruct <-
  readr::read_delim('data/ss_df_wide.csv', delim=',', col_types=cols()) %>%
  select(aa_position, gsp1_aa, consensus_ss)

gsp1_seq <- select(secstruct, gsp1_aa, aa_position) %>% unique() %>% pull(gsp1_aa) %>% as.character()

interfaces <-
  read_delim('data/SASA_interfaces.txt', delim='\t', col_types=cols()) %>%
  filter(partner %in% c('SRM1', 'RNA1')) %>%
  select(-interface) %>%
  unique() %>%
  right_join(data.frame('yeastresnum' = seq(from = 1, to = 220)),
             by = 'yeastresnum') %>%
  complete(partner, yeastresnum) %>%
  filter(!is.na(partner)) %>%
  mutate(in_interface = ifelse((deltarASA < 0.25 | is.na(deltarASA)), F, T)) %>%
  select(-deltarASA)

annotated_sites <- data.frame(
  'site' = c(rep('P-loop', 7), rep('Switch I', 7), rep('Switch II', 9)),
  'position' = c(seq(17,23), seq(39,45), seq(67,75)),
  stringsAsFactors = F)

empiric <-
  readr::read_delim('data/empiric_data_cleaned.txt', delim='\t', col_types=cols()) %>%
  filter(dataset == '15gen') %>%
  left_join(secstruct, by = c('position' = 'aa_position', 'aa_from' = 'gsp1_aa')) %>%
  left_join(interfaces, by = c('position' = 'yeastresnum')) %>%
  mutate(aa_from = factor(aa_from, levels = amino_acid_order),
         aa_to = factor(aa_to, levels = amino_acid_order),
         dataset = factor(dataset, levels = c('6gen', '15gen'))) %>%
  select(mutant, aa_from, position, aa_to, fitness_score, consensus_ss, partner, in_interface) %>%
  left_join(annotated_sites, by = 'position') %>%
  mutate(site = ifelse(!is.na(site), site, 'None'))

mat <-
  empiric %>%
  select(position, aa_to, fitness_score) %>%
  unique() %>%
  pivot_wider(id_cols = aa_to, names_from = position, values_from = fitness_score) %>%
  arrange(desc(aa_to)) %>%
  column_to_rownames('aa_to') %>%
  as.matrix()

mat_noNA <- mat
mat_noNA[is.na(mat_noNA)] <- -12

# cluster the 6gen heatmap
dissim <- get_dist(t(mat_noNA), method = "euclidean")
hc <- hclust(dissim, method = "ward.D2" )

# can rotate using first principle coordinate
pc1 <-
  cmdscale(dissim, eig = T, k = 1)$point %>%
  as_tibble(rownames = 'position') %>%
  arrange(`V1`) %>% pull(position)
hc <- rotate(hc, pc1)

order <- hc$order
pc1
# plot the heatmap
col_fn <- colorRamp2(seq(-4, 4, by=2), RColorBrewer::brewer.pal(5, 'RdBu'))
hm <-
  Heatmap(mat, name='EMPIRIC score\nordered by PC1', col=col_fn, show_heatmap_legend = T, na='black',
          row_title = 'Gsp1 residue', cluster_rows = F, row_names_side = 'left',
          row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_title = 'Gsp1 position',
          # cluster_columns = as.dendrogram(hc),
          # column_order = order,
          column_order = pc1,
          column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          column_names_side = 'top',
          heatmap_legend_param = list('legend_direction' = 'horizontal')
  )

hm_unordered <-
  Heatmap(mat, name='EMPIRIC score', col=col_fn, show_heatmap_legend = T, na='black',
          row_title = 'Gsp1 residue', cluster_rows = F, row_names_side = 'left',
          row_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          row_names_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_title = 'Gsp1 position',
          cluster_columns = F,
          column_title_gp = gpar(fontsize = 8, fontfamily='Helvetica'),
          column_names_gp = gpar(fontsize = 6, fontfamily='Helvetica'),
          column_names_side = 'top',
          heatmap_legend_param = list('legend_direction' = 'horizontal')
  )

GAPGEF_annot = HeatmapAnnotation(
  `GAPGEF` = empiric %>%
    select(position, partner, in_interface) %>%
    unique() %>%
    pivot_wider(names_from = partner, values_from = in_interface) %>%
    column_to_rownames('position') %>%
    as.matrix(),
  col = list(`GAPGEF` = c('TRUE' = 'black', 'FALSE' = 'white')),
  simple_anno_size = unit(3, 'mm'),
  annotation_name_side = 'left',
  border = T, show_legend = F,
  annotation_name_gp = gpar(fontsize = 8)
)

site_annot = HeatmapAnnotation(
  `site` = empiric %>%
    select(position, site) %>%
    unique() %>%
    column_to_rownames('position') %>%
    as.matrix(),
  col = list(`site` = c('P-loop' = 'red', 'Switch I' = 'green', 'Switch II' = 'gold', 'None' = 'white')),
  simple_anno_size = unit(3, 'mm'),
  annotation_name_side = 'left',
  border = T,
  annotation_name_gp = gpar(fontsize = 8)
)


pdf('plots/GAPGEF_interface_labeled_core_only.pdf', width = 16, height = 6)
draw(hm %v% GAPGEF_annot %v% site_annot, heatmap_legend_side = "bottom")
draw(hm_unordered %v% GAPGEF_annot %v% site_annot, heatmap_legend_side = "bottom")
dev.off()




##### ORDER BY GAP GEF kinetics

# prepare kinetics data for plotting as ratio of relative GAP, GEF efficiencies
GAP_kinetics <- read_tsv('data/GAP_kinetics_MichaelisMenten_parameters.txt', col_types = cols()) %>%
  select(mutant, 'kcat_Km' = mean_kcat_Km, 'kcat' = mean_kcat, 'Km' = mean_Km) %>% unique()
GEF_kinetics <- read_tsv('data/GEF_kinetics_MichaelisMenten_parameters.txt', col_types = cols())

WT_GEF <- filter(GEF_kinetics, mutant == 'WT')
WT_GAP <- filter(GAP_kinetics, mutant == 'WT')

GEF_kinetics <-
  GEF_kinetics %>%
  mutate('rel_GEF_kcat_Km' = kcat_Km/WT_GEF$kcat_Km) %>%
  mutate('rel_GEF_kcat' = kcat/WT_GEF$kcat) %>%
  mutate('rel_GEF_Km' = Km/WT_GEF$Km) %>%
  select(mutant, rel_GEF_kcat_Km, rel_GEF_kcat, rel_GEF_Km)

GAP_kinetics <-
  GAP_kinetics %>%
  mutate('rel_GAP_kcat_Km' = kcat_Km/WT_GAP$kcat_Km) %>%
  mutate('rel_GAP_kcat' = kcat/WT_GAP$kcat) %>%
  mutate('rel_GAP_Km' = Km/WT_GAP$Km) %>%
  select(mutant, rel_GAP_kcat_Km, rel_GAP_kcat, rel_GAP_Km)

# EMAP order, taken from clustering in Fig4B of Gsp1 EMAP/APMS manuscript
EMAP_corr_order <- c('D79S', 'T34Q', 'T34E', 'K101R', 'T34G', 'D79A',
                     'T34A', 'Q147E', 'R108I', 'R108L', 'G80A', 'Y157A',
                     'H141E', 'H141R', 'R108Y', 'R108Q', 'R108G', 'Y148I',
                     'H141I', 'R112A', 'R112S', 'R78K')

kinetics <-
  inner_join(GEF_kinetics, GAP_kinetics, by = 'mutant') %>%
  mutate('GAP/GEF' = rel_GAP_kcat_Km/rel_GEF_kcat_Km) %>%
  right_join(data.frame('mutant' = EMAP_corr_order, stringsAsFactors = F)) %>%
  mutate(mutant = factor(mutant, levels = EMAP_corr_order),
         logGAPGEF = log(`GAP/GEF`))

# compute rank order corr between EMPIRIC score and log(GAP/GEF) efficiencies
vec_empiric <-
  empiric %>%
  filter(mutant %in% EMAP_corr_order) %>%
  mutate(mutant = factor(mutant, levels = EMAP_corr_order)) %>%
  arrange(mutant) %>%
  select(mutant, fitness_score) %>%
  unique() %>%
  pull(fitness_score)

vec_GAPGEF <-
  kinetics %>%
  arrange(mutant) %>%
  pull(logGAPGEF)

rank_order_test <- cor.test(vec_empiric, vec_GAPGEF,  method = "spearman")
# pearson_test <- cor.test(vec_empiric, vec_GAPGEF) # r = 0.45, p-val = 0.08


# barplots for log ratios and individual GAP, GEF rel. efficiencies

bp1 <-
  empiric %>%
  filter(mutant %in% EMAP_corr_order) %>%
  mutate(mutant = factor(mutant, levels = EMAP_corr_order)) %>%
  select(mutant, fitness_score) %>%
  unique() %>%
  ggplot(aes(x=fct_rev(mutant), y=fitness_score)) +
  geom_bar(stat='identity', color = 'black', size=0.5) +
  geom_hline(yintercept = -1, linetype='dashed') +
  geom_hline(yintercept = 1, linetype='dashed') +
  xlab('EMPIRIC fitness score') + ylab('Gsp1 mutant') +
  coord_flip()

bp2 <-
  kinetics %>%
  ggplot(aes(x = fct_rev(mutant), y = logGAPGEF)) +
  geom_bar(stat='identity', color = 'black', size=0.5) +
  geom_hline(yintercept = 0) +
  xlab('log(GAP eff / GEF eff)') +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

bp3 <-
  kinetics %>%
  ggplot(aes(x = fct_rev(mutant), y = rel_GAP_kcat_Km)) +
  geom_bar(stat='identity', color = 'black', size=0.5) +
  geom_hline(yintercept = 0) +
  xlab('GAP eff (kcat/Km)') +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

bp4 <-
  kinetics %>%
  ggplot(aes(x = fct_rev(mutant), y = rel_GEF_kcat_Km)) +
  geom_bar(stat='identity', color = 'black', size=0.5) +
  geom_hline(yintercept = 0) +
  xlab('GEF eff (kcat/Km)') +
  coord_flip() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p <- cowplot::plot_grid(bp1, bp2, bp3, bp4, nrow = 1, rel_widths = c(2,2,1,1))

title_text <-
  paste('EMPIRIC score and GAP, GEF kinetic values are correlated',
        paste0("Spearman's rho = ", rank_order_test$estimate %>% round(2)),
        paste0("p-value = ", rank_order_test$p.value %>% round(2)),
        sep = ', ')
title <-
  cowplot::ggdraw() +
  cowplot::draw_label(title_text, x = 0, hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))

cowplot::plot_grid(title, p, ncol = 1, rel_heights = c(0.1, 1))
ggsave('plots/GAP_GEF_kinetics_barplots.png', width = 12, height = 6)
dev.off()




##### Order of all EMAP mutants
# now order by columns in Ext Fig 3
raw_EMAP_order <-
  c('H141E','Y157A','D79A','D79S','T34Q','T34E','K101R','T34G','T34A','R108L',
    'R108I','Q147E','H141I','G80A','H141R','R112A','R112S','R108Y','R108Q',
    'Y148I','R108G','R108A','R78K','K129E','N84Y','A180T','T34D','K129F',
    'T34Y','K129T','T139R','N102K','K143H','R108D','N102M','E115I','Q147L',
    'E115A','F58A','K129I','T137G','T34S','F58L','T139A','K169I','K132H',
    'N105L','N105V','T34L','K143W','H141V','K143Y','R108S','K154M','N102I')

empiric %>%
  filter(mutant == 'A180T')

empiric %>%
  select(mutant, fitness_score) %>%
  unique() %>%
  filter(mutant %in% raw_EMAP_order) %>%
  mutate(mutant = factor(mutant, levels = raw_EMAP_order)) %>%
  ggplot(aes(x=mutant, y=fitness_score)) +
  geom_bar(stat='identity') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('plots/scores_ordered_by_raw_EMAP_clustering.png', width = 12, height = 3)
