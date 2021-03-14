import pysam
import sys
import pandas as pd

def is_split_read(s):
    try:
        s.get_tag("SA")
    except KeyError:
        return False
    return True


def get_read_map(bam_path, contig_of_interest='NC_007605'):
    coordinates = pd.DataFrame()

    infile = pysam.AlignmentFile(bam_path, "rb")
    for s in infile.fetch(contig_of_interest):
    
        if s.mapping_quality<30 or s.is_duplicate or \
            s.is_unmapped or s.mate_is_unmapped: 
            continue            
            
        if is_split_read(s):
            sa = s.get_tag("SA").split(",")
            chrom1, pos1 = sa[0], int(sa[1])
        elif not s.is_proper_pair:               
            chrom1, pos1 = s.next_reference_name, s.next_reference_start            
        else: # ignore normal mapping reads
            continue

        if chrom1 != contig_of_interest:
                        
            coordinates = coordinates.append({'readpair': s.query_name, 'chr1': chrom1, 'pos1': pos1, 
                                              'chr2': s.reference_name, 'pos2': s.reference_start}, 
                                             ignore_index=True)    
            
    return coordinates[['readpair', 'chr1', 'pos1', 'chr2', 'pos2']].astype({'pos1':int, 'pos2':int})


def get_clusters_df(df, chr_col, pos_col, cluster_colname, cluster_name_prefix,
                    MAX_DIST=1000, MIN_TAGS=6, drop_no_clust=False):
    """
    Cluster by chr_col and pos_col. Reads are clustered together if there is at least MIN_TAGS reads pairwise separated by MAX_DIST.
    So, with default params, reads in positions 1000,2000,3000,4000,5000,and 6000, will be enough to compose a cluster.
    Cluster ID, composed of cluster_name_prefix and a sequential integer, is placed into cluster_colname for the reads in the cluster. 
    Not clustered reads will have the value in cluster_colname unchanged.
    """
    res=df.copy().sort_values([chr_col, pos_col])
    ccnt = 1
    cstart=0
    clust_num = 1
    for i in range(1,len(res)):
        if res[pos_col].iloc[i] - res[pos_col].iloc[i-1] <= MAX_DIST:
            ccnt += 1
        else:
            if ccnt >= MIN_TAGS:
                res[cluster_colname].iloc[range(cstart, i)] = (cluster_name_prefix + str(clust_num))
                clust_num += 1
            ccnt = 1
            cstart = i
    # print last cluster
    if ccnt > MIN_TAGS:
        res[cluster_colname].iloc[range(cstart, i)] = (cluster_name_prefix + str(clust_num))
    
    if drop_no_clust:
        return(res[res[cluster_colname]!=""])
    else:
        return res

    
if __name__ == '__main__':
    
    coordinates = get_read_map(sys.argv[1], 'NC_007605')
    coordinates['integr_hotspot'] = ''
    coordinates['EBV_cluster'] = ''

    max_distance = 1000
    min_tags = 6

    final_df = pd.DataFrame()
    for chr_ in coordinates['chr1'].unique():

        chr_df = get_clusters_df(coordinates[coordinates['chr1']==chr_], 
                                 chr_col='chr1', pos_col='pos1',
                                 cluster_colname='integr_hotspot',
                                 cluster_name_prefix=chr_+'_hs',
                                 MAX_DIST = max_distance, 
                                 MIN_TAGS = min_tags)
        if len(chr_df)<=0: 
            continue

        # add readpairs without hotspot
        final_df = pd.concat([final_df, chr_df[chr_df['integr_hotspot'] == '']])

        # cluster ebv_reads within each integration hotspot
        for hs in [e for e in chr_df['integr_hotspot'].unique() if e]:
            chr_hs_df = get_clusters_df(chr_df[chr_df['integr_hotspot'] == hs], 
                                 chr_col='chr2', pos_col='pos2',
                                 cluster_colname='EBV_cluster',
                                 cluster_name_prefix=hs+'_ebv_cluster',
                                 MAX_DIST = max_distance, 
                                 MIN_TAGS = min_tags)
            final_df = pd.concat([final_df, chr_hs_df])

    final_df = final_df[['readpair', 'integr_hotspot', 'chr1', 'pos1', 'chr2', 'pos2', 'EBV_cluster']].sort_values(['chr1','pos1','pos2'])
    
    final_df.to_csv(sys.stdout, sep='\t', index=False)
