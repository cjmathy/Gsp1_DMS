#!/usr/bin/env RScript

# run from command line using 'Rscript empiric_preprocessing.R' from
# the scripts directory of the project

library(tidyverse)

# Store data file names in a named list, making reading into a df easier
datafiles <- list('6gen' = 'data/GSP1_6gen_fitness_20181105.csv',
                  '15gen' = 'data/GSP1_15gen_fitness_20181105.csv')

# set ordering for amino acids
amino_acid_order <- c('*','K','R','H','D','E','N',
                      'Q','S','T','P','G','C','A',
                      'V','I','L','M','F','Y','W')

# read in the WT sequence for Gsp1 and make it a dataframe for merging
# FASTA file downlaoded from UNIPROT on 2020-01-20
# https://www.uniprot.org/uniprot/P32835.fasta
sequence <-
  read.delim('data/gsp1.fasta', header = T, stringsAsFactors = F)[[1]] %>%
  paste(collapse='') %>%
  strsplit(split='')
sequence[[1]][220] = '*'  # need to add STOP codon code at position 220

seq_df <-
  data.frame(sequence = sequence,
             'position' = seq_len(length(sequence[[1]])),
             stringsAsFactors = F, fix.empty.names = F)
names(seq_df) <- c('sequence', 'position')


# read in the files, correctly name the columns, and merge the WT sequence
empiric <-
  map_dfr(datafiles, read_csv,
          .id = 'dataset', na = '--',
          col_names = c('mutation','fitness_score'),
          col_types = cols(mutation = col_character(),
                           fitness_score = col_double())) %>%
  separate(mutation, c('mutated_to','position'), sep = 1) %>%
  mutate('aa_to' = factor(mutated_to, levels = amino_acid_order),
         'position' = as.numeric(position)) %>%
  left_join(seq_df, by='position') %>%
  mutate(mutant = paste0(sequence, position, aa_to),
         'aa_from' = factor(sequence, levels = amino_acid_order)) %>%
  select(mutant, aa_from, position, aa_to, dataset, fitness_score)

write_delim(empiric, path='data/empiric_data_cleaned.txt', delim='\t')

empiric %>%
  filter(dataset == '15gen') %>%
  select(mutant, aa_from, position, aa_to, fitness_score) %>%
  arrange(desc(fitness_score), position) %>%
  write_csv(path='data/15gen_mutation_scores_ordered.csv')
