# 2020-11-13 Aggregation Test

- On 2020-11-10, Tina purified F56V and F58Y (frozen pellets from 1L of culture).
- I concentrated them and froze 3 x 70 uL aliquots each, at 104 uM (F58V) and 300 uM (F58Y) 
- On 2020-11-13, Tina loaded the mutants to test if they aggregate during GTP Loading. I ran HPLC and CD on them

## HPLC

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import glob

expt_dir = '../data/HPLC/20201112_GSP1_F28V_F28Y_S75 2020-11-12 18-26-05/'
datafiles = glob.glob(expt_dir+'*/*280NM.CSV',recursive=True)

dfs = []

for f in datafiles:
  df = pd.read_csv(f, encoding='utf-16-le', sep='\t', names=('retention (mL)','nm280 (mAU)'))
  df['sample'] = f.split('/')[-1].split('.')[0]
  dfs.append(df)

df = pd.concat(dfs)

sns.relplot(data=df, x='retention (mL)', y='nm280 (mAU)', col='sample', kind='line',
            col_order=['BLANK1_SEC75_280NM','F28V_SEC75_280NM','BLANK2_SEC75_280NM','F28Y_SEC75_280NM'])

