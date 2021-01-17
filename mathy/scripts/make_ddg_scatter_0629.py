import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

%matplotlib inline
sns.set(style="whitegrid", palette="pastel", color_codes=True, font='Helvetica')

# Set font sizes
SMALL_SIZE = 12 #6
MEDIUM_SIZE = 13 #7
BIG_SIZE = 14 #8

plt.rcParams.update({
    'legend.loc': 'center left',
    'font.family': 'Helvetica',
    'font.size': BIG_SIZE,         # controls default text sizes
    'axes.titlesize': SMALL_SIZE,  # fontsize of the axes title
    'axes.labelsize': MEDIUM_SIZE, # fontsize of the x and y labels
    'xtick.labelsize': SMALL_SIZE, # fontsize of the tick labels
    'ytick.labelsize': SMALL_SIZE, # fontsize of the tick labels
    'legend.fontsize': SMALL_SIZE, # legend fontsize
    'figure.titlesize': BIG_SIZE,   # fontsize of the figure title
    'figure.figsize': (10,10),
    'figure.dpi': 300
})

ddg_dir = '../ddg_monomer_ala_scan_20200629/'

# read in ddg data
ddg_df = pd.read_csv(ddg_dir+'ddg_results.txt', delim_whitespace=True)

# read in fitness data
dms_df = pd.read_csv('../../Data/15gen_binned_0607.csv')
dms_df.bin = dms_df.bin.mask(pd.isnull(dms_df.bin), other='drop-outs')
bin_order = ['drop-outs', 'toxic', 'STOP-like', 'intermediate', 'neutral', 'beneficial']


# fix numbering, because "description" lists mutations using re-numbered
# residues for 3m1i with the tail chopped off, rosetta numbered residues
# 1-174 correspond to Gsp1 residues 10-183

ddg_df['aa_from'] = ddg_df['description'].str[0]
ddg_df['aa_to'] = ddg_df['description'].str[-1]
ddg_df['rosetta_resnum'] = ddg_df['description'].str[1:-1].astype(int)
ddg_df['position'] = ddg_df['rosetta_resnum']+9


# merge the datasets
# this will only have 173 out of 174 calculated ddGs, because of Q71L
# mutation in the structure (so we computed L71A but the DMS scan did Q71A)

df = pd.merge(dms_df, ddg_df, how='inner', on=['aa_to','position','aa_from'])
df.rename(columns={'score': 'fitness', 'total': 'ddG'}, inplace=True)
df.drop(['description','rosetta_resnum'], axis=1, inplace=True)


# make NaN's in fitness equal the minimum value, so we see them on the plot
df.fitness = df.fitness.fillna(min(df.fitness))
df.bin = df.bin.apply(lambda row: 'drop-outs' if pd.isnull(row) else row)


# plot
fig, ax = plt.subplots()
sns.scatterplot(data=df, x='fitness', y='ddG', hue='bin',
                hue_order = bin_order, palette=['black']+sns.color_palette('RdBu', n_colors=7)[:5],
                edgecolor='black', linewidth=0.5)
plt.title('Alanine scan, residues 10-183 (no tail), PDB=3m1iA')
plt.xlabel('Fitness Score')
plt.ylabel('ddG (REU - Rosetta Energy Units)')
plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
fig.set_size_inches(10, 8)
fig.tight_layout(rect=[0, 0, 1, 1])
fig.savefig('../Mathy_plots/alascan_ddG_2020-06-29.png', dpi=300)

