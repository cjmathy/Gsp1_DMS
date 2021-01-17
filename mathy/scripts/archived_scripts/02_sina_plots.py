#!/usr/bin/env python3

# This script defines groups of positions based on known structural
# annotation, then plots sina/violin plots of distributions to directly
# compare the selection coefficients of mutations at positions belonging
# to different groups.

# cjmathy@gmail.com 2020-04-27

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sinaplot import sinaplot

WT_seq = 'MSAPAANGEVPTFKLVLVGDGGTGKTTFVKRHLTGEFEKKYIATIGVEVHPLSFYTNFGEIKFDVWDTAGQEKFGGLRDGYYINAQCAIIMFDVTSRITYKNVPNWHRDLVRVCENIPIVLCGNKVDVKERKVKAKTITFHRKKNLQYYDISAKSNYNFEKPFLWLARKLAGNPQLEFVASPALAPPEVQVDEQLMQQYQQEMEQATALPLPDEDDADL*'

df = pd.read_csv('../Data/15gen_binned.csv')

# define functional groups based on residue positions

func_groups = {
    'P-loop': list(range(17, 23+1)),
    'Switch I':  list(range(39, 45+1)),
    'Switch II': list(range(67, 75+1)),
    'C-terminal linker': list(range(181,193+1)),
    'C-terminal helix': list(range(194,208+1)),
    'Effector lobe': list(range(1, 96+1)),
    'Allosteric lobe': list(range(97,180+1)),
    'C-terminal extension': list(range(181, 220+1))
}

func_groups_df = (
    pd.DataFrame.from_dict(data=func_groups, orient='index')
    .transpose()
    .melt(var_name='region', value_name='position')
    .dropna()
    .assign(position=lambda x: x['position'].astype('int'))
    .reindex(columns=['position','region'])
)

df = df.merge(func_groups_df, on='position')


sns.set(style="whitegrid", palette="pastel", color_codes=True)


df1 = df[df.region.isin(['P-loop','Switch I','Switch II', 'C-terminal linker',
                         'C-terminal helix'])]
df1 = df1[df1.bin.isin(['beneficial'])]


sns.catplot(x='region', y='score', kind='violin',inner=None,data=df1)



# Draw a nested violinplot and split the violins for easier comparison
# sns.violinplot(x='bin', y='score', hue='region',
 #              split=True, inner="quart",
 #              palette={'Effector lobe':'c', 'Allosteric lobe': 'g'},
 #              data=df[df.region.isin(['Effector lobe','Allosteric lobe'])])
 #sns.despine(left=True)
