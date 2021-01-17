#!/usr/bin/env python
# coding: utf-8

import pandas as pd
# split conservation and SASA information

in_file = '~/gdrive/gsp1_dms/mathy/data/data_from_tina/consurf_and_mds_and_sasa.txt'
out_file_burial = '~/gdrive/gsp1_dms/mathy/data/burial.csv' 
out_file_conserv = '~/gdrive/gsp1_dms/mathy/data/conservation.csv' 

df = pd.read_csv(in_file, delimiter='\t')

(df[['position','aa_to','consurf','consurf_score','exists_in_msa']]
 .drop_duplicates()
 .to_csv(out_file_conserv, index=False)
)

(df[['position','region3','sasa_group_category']]
 .drop_duplicates()
 .to_csv(out_file_burial, index=False)
)
