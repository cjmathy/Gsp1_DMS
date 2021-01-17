import pandas as pd

index = pd.read_csv('../Data/gene_sets/gene_sets_with_ORF.txt', sep='\t')
clustering_gene_sets = pd.read_csv('../Data/gene_sets/corr_clustering_heatmap_cluster_gene_sets.csv')
df = pd.merge(clustering_gene_sets, index, left_on='strain',right_on='query')
df.rename(columns={'cluster_x':'set_clustered_hmap',
                   'cluster_y':'set_clustered_tina',
                   'gene_set':'set_handmade_tina',
                   'query_ORF':'ORF',
                   'gene_name':'name'},
          inplace=True)
df = df[['strain', 'name', 'ORF', 'set_clustered_hmap','set_clustered_tina',
         'set_handmade_tina','Description']]

df.to_csv('../Data/gene_sets/heatmap_gene_sets_2020_0505.csv')

# alternatively, just just tina's "gene_sets_final.csv"

# from this file, we want the strain name (query), the gene_name, the gene_set,
# and the cluster
final_gene_sets = pd.read_csv('../Data/gene_sets/gene_sets_final.txt', sep='\t')

# from this file, we want the cluster descriptions (i.e. pos_with_GAP) and the descriptions
tina_gene_sets = (pd.read_csv('../Data/gene_sets/gene_sets_with_ORF.txt', sep='\t')
                  [['query','cluster','Description','query_ORF']])

# from this file we want the ordering of the heatmap columns
clustering_gene_sets = (pd.read_csv('../Data/gene_sets/corr_clustering_heatmap_cluster_gene_sets.csv')
                        .reset_index(level=0)
                        .rename(columns={'index':'strain_order'})
                        .assign(cluster=lambda x: x['cluster'].astype('object'))
                        [['strain_order','strain']]
                       )

# now merge tina_gene_sets and clustering_gene_sets into final_gene_sets

df = (pd.merge(final_gene_sets, tina_gene_sets, on='query', suffixes=('','_y'))
      .drop_duplicates()
      .rename(columns={'cluster_y':'set_info',
                       'query':'strain',
                       'gene_name':'name',
                       'query_ORF':'ORF'}
             )
     )

# now merge in the cluster ordering
df = (pd.merge(df, clustering_gene_sets, how='left', on='strain')
      [['strain','name','ORF','gene_set',
        'set_info','cluster','strain_order','Description']]
      .sort_values(by=['cluster','strain_order','gene_set','name'])
      .assign(strain_order=lambda x: x['strain_order'].astype('Int64')+1)
      .reset_index(drop=True)
     )

df.to_csv('../Data/gene_sets/complete_gene_sets_2020_0505.csv')
