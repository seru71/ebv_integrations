# ebv_integrations

Script indentifying pairs of reads that map to EBV genome (NC_007605) and one of human chromosomes, and clustering them into integration hotspots (on chromosomes) and donor-clusters (on EBV genome).

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
