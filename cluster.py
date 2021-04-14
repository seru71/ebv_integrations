import pandas as pd
import sys

infile = sys.argv[1]
outfile = infile[:-len('clustered_reads.tsv')] + 'clusters.tsv'
print('Reading', infile)
t = pd.read_csv(infile, sep='\t')

# get fragments clustered on both ends
t = t[t['integr_hotspot'].notna()]
t = t[t['EBV_cluster'].notna()].iloc[:,1:7]


# get hotspot coords
t1 = t.groupby(['integr_hotspot']).agg(['min','max']).reset_index()
t1 = pd.DataFrame({'hotspot': t1.iloc[:,0],
                   'chr':   t1.loc[:,('chr1', 'min')],
                   'start': t1.loc[:,('pos1', 'min')],
                   'end':   t1.loc[:,('pos1', 'max')],
                  })

# get EBV coords
t2 = t.groupby(['EBV_cluster']).agg(['min','max','count']).reset_index()
t2 = pd.DataFrame({'hotspot': t2.iloc[:,1],
                   'ebv_cluster': t2.iloc[:,0],
                   'ebv_chr':   t2.loc[:,('chr2', 'min')],
                   'ebv_start': t2.loc[:,('pos2', 'min')],
                   'ebv_end':   t2.loc[:,('pos2', 'max')],
                   'fragments': t2.loc[:,('pos2', 'count')],
                  })
                  
# merge
result = pd.merge(left=t1, right=t2, left_on='hotspot', right_on='hotspot')

# save
print('Writing', outfile)
result.to_csv(outfile, sep='\t', index=False)
