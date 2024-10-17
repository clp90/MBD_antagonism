# MBD_antagonism
Scripts and code for MBD5/6 and MBD7 antagonism paper (2024) - link **TBD**

General Dependencies:
- all analyses for this project used genome annotations that we previously updated based on pollen RNA-seq data, the GTF files are available here: https://github.com/clp90/mbd56_pollen
- all sequencing data have been deposited to GEO under SuperSeries GSE275832 

## snmCT-seq data analysis
### snmCT_seq_map.sh
This script handles all initial preprocessing, mapping, and basic postprocessing for one 384 well plate of snmCT-seq data. Dependencies are listed below. Uses `bismark` for WGBS mapping and `STAR` for RNA-seq mapping, so genome indexes are required for both of these. Genome annotations in GTF format are also required. Requires a list of barcodes provided in `snmCAT_well_barcodes.txt` in the helper_scripts/ subdir. See header of this script for examples, or run the script with the -h option to get usage info. Run with -0 to test that all dependencies are installed.  

Example usage:
```
snmCT_seq_map.sh -1 forward_reads.fastq.gz [-2 reverse_reads.fastq.gz] -b bismark_index -t STAR_index -o output_folder -i helper_scripts/snmCAT_well_barcodes.txt -g annotations.gtf [options]
```
- Must be installed on your `$PATH`:
	- fastqc (from Babraham Institute, tested on v0.11.8)
	- trim_galore (from Babraham Institute, tested on v0.4.1)
	- STAR (by Dobin et al. 2013, tested on v2.7.9a)
	- bismark (by Babraham Bioinformatics, v.0.20.1)
	- bowtie2 (by Ben Langmead, v.2.3.4.3)
	- python (v3.9.6)
	- HTseq (Anders et al. 2014 Bioinformatics, tested on v2.0.2) - to install: pip install --user HTseq
	- cutadapt (v.3.4 by Marcel Martin) - to install: pip install --user cutadapt
	- trim_galore (from Babraham Institute, tested on v0.6.7)
	- samtools (by Heng Li, tested on v1.11)
	- bedGraphToBigWig (part of UCSC tools, tested on v4)
	- htseq-count (part of HTSeq python package, tested on v0.13.5)
	- bedtools (v.2.30.0 by quinlanlab.org and others)
	- java (tested on v.1.7.0_181, OpenJDK Runtime Environment (IcedTea 2.6.14) (7u181-2.6.14-0ubuntu0.2), OpenJDK 64-Bit Server VM (build 24.181-b01, mixed mode))
- Must be in helper_scripts/ subdirectory, or set path to this folder using -s option:
	- deduplicate_bismark.pl (by Babraham Bioinformatics, downloaded Dec. 18, 2018)
	- bismark_summary_to_bed.py (by Colette Picard, v.1.0, 10/26/2014)
	- steve_filter.py (v.1.0 by Colette Picard, 02/28/2020)
- Other requirements:
	- Markduplicates from picard_tools v.2.25.0 : set path to picard_tools jar using -P (default assumes it's saved to a variable named PICARD)

### snmCAT_seq_QC_plots.R
This is a simple script that makes a figure to visualize the final QC status of each well of a 384 well plate processed using `snmCT_seq_map.sh`. This can help find issues with sorting (e.g. when all wells after a certain point are empty/failed QC) or other library prep issues.  
Example usage:
```
snmCAT_seq_QC_plots.R snmCT_seq_map_outputdir/summaries/all_statuses.txt outfilename.pdf --mode 1
```
For the moment 'mode 1' is the only supported option, this may be extended in future to make other QC plots.
- Must be installed on your `$PATH`:
  - R (v.4.1.1)
    - R packages `argparse`, `ggplot2`, `reshape2`

## Pseudobulking snmCT-seq WGBS data
### pseudobulk_WGBS.sh
This script pools WGBS data across all wells of a 384-well plate according to a list of nuclei cluster assignments. Takes the output of snmCT_seq_map.sh and a list of cluster assignments as inputs. See header of this script for examples, or run the script with the -h option to get usage info. Run with -0 to test that all dependencies are installed.  
Output is a BED-like file with per-position data, seven fields: chromosome, start, end, # unmethylated reads, # methylated reads, percent methylation, number of nuclei that had coverage over this site.  
Example usage:
```
pseudobulk_WGBS.sh -d snmCTseq_dir/WGBS_per_position/ -c cluster_assn.txt -o output_folder -f genome.chrom.sizes [options]
```
Note that the chrom.sizes file will be in the snmCT_seq_map.sh output directory.  
Required on your `$PATH`:
- bedtools (v.2.30.0 by quinlanlab.org and others)
- bedGraphToBigWig (part of UCSC tools, tested on v4)

### pool_WGBS_bed.sh
This is a small helper script that will pool the output of multiple runs of `pseudobulk_WGBS.sh`, allowing you to pool data across multiple plates.  
Output is a BED-like file with per-position data, seven fields: chromosome, start, end, # unmethylated reads, # methylated reads, percent methylation, number of nuclei that had coverage over this site, number of input files that had coverage over this site.
Example usage:
```
pool_WGBS_bed.sh output_prefix inputfile1.bed inputfile2.bed ...
```
Required on your `$PATH`:
- bedtools (v.2.30.0 by quinlanlab.org and others)

## Calculating weighted average methylation over genomic regions
### sumByFeature.sh
This is a small script that calculates the weighted mean methylation over a region (sum of methylated counts / sum of methylated and unmethylated counts, across all sites overlapping the region). Input is expected to be a BED-like file of per-position methylation data containing six fields: chromosome, start, end, # unmethylated reads, # methylated reads, percent methylation, and a BED file of intervals that you want to calculate average methylation over. If your input files are both sorted by position (sort -k1,1 -k2n,2) then you can add the -S option to make this run significantly faster.
Example usage:
```
sumByFeature.sh -i methylation_data.bed -r regions.bed -o outputfile.bed [options]
```
Required on your `$PATH`:
- bedtools (v.2.30.0 by quinlanlab.org and others)
- python v.3, tested on v.3.9.6
Other requirements (must be in same location as this script):
- summarizeMethylation.py by CLP




