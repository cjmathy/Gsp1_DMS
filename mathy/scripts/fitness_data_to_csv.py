#!/usr/bin/env python3

# This script reads in the pickled files with scores and bins, and then
# preprocesses the data into a long-form data table, saved as a csv.
# This csv is used for the further analysis

# chris.mathy@ucsf.edu
# last updated 2020-09-29

import numpy as np
import pandas as pd
import pickle as pkl

# amino acid integer code
transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5,
             'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15,
             'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: '*', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M',
                 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T',
                 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}

# bin float code
bins_dict = {0.5: 'beneficial',
             0: 'neutral',
             -0.5: 'intermediate',
             -1: 'STOP-like',
             -2: 'toxic'}

# Prepare a dataframe with the WT sequence for adding to the
# mutant data
# For dataframes where each row is a mutation, the WT residues
# is found in the column 'aa_from'

WT_seq = 'MSAPAANGEVPTFKLVLVGD' \
         'GGTGKTTFVKRHLTGEFEKK' \
         'YIATIGVEVHPLSFYTNFGE' \
         'IKFDVWDTAGQEKFGGLRDG' \
         'YYINAQCAIIMFDVTSRITY' \
         'KNVPNWHRDLVRVCENIPIV' \
         'LCGNKVDVKERKVKAKTITF' \
         'HRKKNLQYYDISAKSNYNFE' \
         'KPFLWLARKLAGNPQLEFVA' \
         'SPALAPPEVQVDEQLMQQYQ' \
         'QEMEQATALPLPDEDDADL*'

WT_seq_df = (
    pd.DataFrame({'aa_from': [res for res in WT_seq]})
    .rename_axis('position')
    .reset_index()
)

### 15 generation dataset to CSV ### 

# read in the scores and binning labels from pkl files
score_file = '../../Data/15Gen_presort_sorted.pkl'
bins_file = '../../Data/15Gen_presort_sorted_bins.pkl'
score_mat_15 = pkl.load(open(score_file, 'rb'), encoding='latin1')
bins_mat_15 = pkl.load(open(bins_file, 'rb'), encoding='latin1')

# Create dataframes in long format
scores = (
    pd.DataFrame(data=score_mat_15)
    .rename_axis('aa_to')
    .reset_index()
    .pipe(pd.melt, var_name='position', value_name = 'score', id_vars='aa_to')
)

bins = (
    pd.DataFrame(data=bins_mat_15)
    .rename_axis('aa_to')
    .reset_index()
    .pipe(pd.melt, var_name='position', value_name = 'bin', id_vars='aa_to')
    .assign(bin=lambda x: x['bin'].astype('float').map(bins_dict))
    .assign(bin=lambda x: x['bin'].apply(lambda b: 'drop-outs' if pd.isnull(b) else b))
)

# merge the dataframes
df = (
    pd.merge(scores, bins, on=['aa_to','position'])
    .assign(aa_to=lambda x: x['aa_to'].map(transform_inv))
    .merge(right=WT_seq_df, on='position')
    .assign(position=lambda x: x['position']+1)
    .assign(mutant=lambda x: x['aa_from']+x['position'].astype(str)+x['aa_to'])
    .reindex(columns=['mutant','aa_from','position','aa_to','score','bin'])
)

# save the 15 generation dataset
df.to_csv('../../Data/15gen_fitness_current.csv', index=False, na_rep='NaN')

### 6 generation dataset to CSV ### 

# read in the scores and binning labels from pkl files
score_file = '../../Data/6gen_data_matrix.pkl'
binning_file = '../../Data/6gen_data_matrix_bins.pkl'
score_mat_6 = pkl.load(open(score_file, 'rb'), encoding='latin1').filled(fill_value=np.nan)
bins_mat_6 = pkl.load(open(binning_file, 'rb'), encoding='latin1').filled(fill_value=np.nan)

# Create dataframes in long format
scores = (
    pd.DataFrame(data=score_mat_6)
    .rename_axis('aa_to')
    .reset_index()
    .pipe(pd.melt, var_name='position', value_name = 'score', id_vars='aa_to')
)

bins = (
    pd.DataFrame(data=bins_mat_6)
    .rename_axis('aa_to')
    .reset_index()
    .pipe(pd.melt, var_name='position', value_name = 'bin', id_vars='aa_to')
    .assign(bin=lambda x: x['bin'].astype('float').map(bins_dict))
    .assign(bin=lambda x: x['bin'].apply(lambda b: 'drop-outs' if pd.isnull(b) else b))
)

# merge the dataframes
df = (
    pd.merge(scores, bins, on=['aa_to','position'])
    .assign(aa_to=lambda x: x['aa_to'].map(transform_inv))
    .merge(right=WT_seq_df, on='position')
    .assign(position=lambda x: x['position']+1)
    .assign(mutant=lambda x: x['aa_from']+x['position'].astype(str)+x['aa_to'])
    .reindex(columns=['mutant','aa_from','position','aa_to','score','bin'])
)

# save the 15 generation dataset
df.to_csv('../../Data/6gen_fitness_current.csv', index=False, na_rep='NaN')