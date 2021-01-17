#!/usr/bin/env python3

# This script reads in the pickled files with scores and bins, and then
# preprocesses the data into a long-form data table, saved as a csv.
# This csv is used for the further analysis

# cjmathy@gmail.com 2020-04-27

import numpy as np
import pandas as pd
import pickle as pkl

transform = {'A': 9, 'C': 8, 'E': 20, 'D': 19, 'G': 10, 'F': 2, 'I': 5,
             'H': 16, 'K': 18, '*': 0, 'M': 6, 'L': 4, 'N': 14, 'Q': 15,
             'P': 11, 'S': 12, 'R': 17, 'T': 13, 'W': 1, 'V': 7, 'Y': 3}
transform_inv = {0: '*', 1: 'W', 2: 'F', 3: 'Y', 4: 'L', 5: 'I', 6: 'M',
                 7: 'V', 8: 'C', 9: 'A', 10: 'G', 11: 'P', 12: 'S', 13: 'T',
                 14: 'N', 15: 'Q', 16: 'H', 17: 'R', 18: 'K', 19: 'D', 20: 'E'}
bins_dict = {0.5: 'beneficial', 0: 'neutral', -0.5: 'intermediate',
             -1: 'STOP-like', -2: 'toxic'}

WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'

score_matrix = pkl.load(open('../Data/15_Gen_with_missing.pkl', 'rb'),
                        encoding='latin1')
bins_matrix = pkl.load(open('../Data/15_Gen_with_missing_bins.pkl', 'rb'),
                       encoding='latin1')

scores = (
    pd.DataFrame(data=score_matrix)
    .rename_axis('aa_to')
    .reset_index()
    .pipe(pd.melt, var_name='position', value_name = 'score', id_vars='aa_to')
)


bins = (
    pd.DataFrame(data=bins_matrix)
    .rename_axis('aa_to')
    .reset_index()
    .pipe(pd.melt, var_name='position', value_name = 'bin', id_vars='aa_to')
    .assign(bin=lambda x: x['bin'].astype('float').map(bins_dict))
)

WT_seq_df = (
    pd.DataFrame({'aa_from': [res for res in WT_seq]})
    .rename_axis('position')
    .reset_index()
)

df = (
    pd.merge(scores, bins, on=['aa_to','position'])
    .assign(aa_to=lambda x: x['aa_to'].map(transform_inv))
    .merge(right=WT_seq_df, on='position')
    .assign(position=lambda x: x['position']+1)
    .assign(mutant=lambda x: x['aa_from']+x['position'].astype(str)+x['aa_to'])
    .reindex(columns=['mutant','aa_from','position','aa_to','score','bin'])
)

df.to_csv('../Data/15gen_binned.csv', index=False, na_rep='NaN')
