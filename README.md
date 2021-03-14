# ebv_integrations

Script indentifying pairs of reads that map to EBV genome (NC_007605) and one of human chromosomes, and clustering them into integration hotspots (on chromosomes) and donor-clusters (on EBV genome).

## Pipeline

1. Select reads mapping to NC_007605 with MQ>=30, and being either:
  
   a) split between NC_007605 and another contig
  
   b) discordant with mate mapping to reference contig other than NC_007605
  
2. Cluster reads based on mapping position to the chromosome (integration hotspots). Reads are clustered if:

   a) distance between consecutive mapping positions is <= 1kb, and
   
   b) at least 6 such reads are found
   
3. Within each integration hotspot, cluster mapping positions on EBV genome to identify "donor clusters". 
   The same criteria as above are used to group reads into clusters.

## Usage

`python3 get_cluster_pairs.py BAM_FILE > result.tsv`

## Result

Table with fragments (read pairs) where one of the mates maps to EBV genome (NC_007605) and the other to other contig in the reference genome.
Table columns:
 - `readpair` - id of the fragment / read-pair
 - `integr_hotspot` - id of the hotspot, incuding name of the chromosome and a consecutive number. If empty the read was not in an integration hotspot
 - `chr1` - contig where one of the mates mapped other than NC_007605
 - `pos1` - position on chr1 where the read mapped
 - `chr2` - contig NC_007605
 - `pos2` - position on NC_007605 where the read mapped 
 - `EBV_cluster` - id of the "donor cluster", incuding name of the linked chromosome and a consecutive number. If empty the read was not included in any cluster on EBV genome

## Notes

Read counts in the clusters should be normalized by the total number of reads in a given sample. 
Data for it can be found in `data/*.bas` files.

## Requirements
 - pysam
 - pandas
