#!/usr/bin/env bash
trap 'pkill -P $$; exit' SIGINT SIGTERM					# kills child processes on exit

# ------------------------------------------------------------------------------------
# v3.1 by Colette Picard
# 03/15/2022
# ------------------------------------------------------------------------------------
	
# Usage:
# scmCT_seq_map.sh [options] -1 forward.fq -g bismark_index -t STAR_index -o outdir

# -------------------------
# Version history:
# v.1.0: initial build - 03/15/2022 by CLP
# -------------------------


# Description printed when "help" option specified:
read -d '' usage <<"EOF"
v1.0 by Colette Picard, 03/15/2022
This script runs a pipeline for the initial processing and aligning of data from
scmCAT-seq, a single-cell multiome RNA-seq + WGBS plate-based method published here:
https://www.cell.com/cell-genomics/fulltext/S2666-979X(22)00027-1

Luo et al. (2022) Single nucleus multi-omics identifies human cortical cell regulatory
genome diversity. Cell Genomics, DOI: https://doi.org/10.1016/j.xgen.2022.100107

OVERVIEW:
(Step 1) Well-level demultiplexing
(Step 2) Read trimming and quality check
(Step 3) Align reads to transcriptome using STAR
(Step 4) Align reads to bisulfite-treated genome using bismark
(Step 5) Check library complexity + remove PCR duplicates
(Step 6) Post-processing of RNA-seq data (read counting, etc.)
(Step 7) Post-processing of WGBS data (methylation extractor, getting average counts over genome tiles)
(Step 8) Final summary and cleanup

-------------------------------
ADDITIONAL CONSIDERATIONS:
- This script highly benefits from having access to multiple threads; where possible the
wells will be analyzed in parallel across all available threads.
- Can also be used for snBS-seq or snRNA-seq alone done using the same plate-based format
- Reads can be paired-end or single-end. Paired-end reads are treated as unpaired, however,
because the library prep method generates a high rate of chimeras (10% or more). Well index
is the first 8bp of read1 if paired-end.
-------------------------------

Expected format for barcode index file (-i) is FASTA, for example:
>A1
^ACGATCAG
>A3
^TCGAGAGT
>A5
^CTAGCTCA

(* the '^' at the beginning of the line indicates that each barcode is the first 8bp)


EXTRA STUFF:
-------------------------------
TBD
-------------------------------

Usage:
snmCAT_seq_map.sh [options] -1 forward.fq -g bismark_index -t STAR_index -o outdir

User-specified options (defaults indicated in brackets):
Required arguments:
	-1 forward : name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok)
	-b bismark_index : bismark-indexed reference genome for alignment (see bismark_genome_preparation)
	-t STAR_index : STAR-indexed transcriptome for alignment
	-g gtf : gene annotation file in GTF format
	-o outdir : name of directory with all output files; will be created if nonexistant
Additional options:
	-2 reverse : name of file containing reverse reads if using paired-end data (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok) []
	-n name : short experiment name used to name output files ["expt"]
	-s path_to_scripts : path to folder containing all required scripts (note - $scriptDir in default is the location of this script) ["$scriptDir/scripts"]
	-S path_to_scratch : path to folder to use for scratch/fast calculations; WARNING: FILES AT THIS LOCATION MAY BE OVERWRITTEN, do not point to important folders (will be in $outdir by default, but you can speed up this script by pointing to a disk w/ faster I/O) []
	-p threads : indicate how many threads are available to this program (note: bismark uses 4 threads by default, recommend using at least 4) [1]
	-i idxfile : (demux) path to file containing well barcodes in fasta format ["$scriptDir/snmCAT_well_barcodes.txt"]
	-c trim : trims t bases from 5' and 3' end of reads before any additional quality- or adapter-based trimming (for read w/ well barcode, additional 8bp is added to this); recommend setting to 10 given adaptase behavior (see trim_galore option -t) [10]
	-q minQ : (trim_galore) min quality to use for trimming off bad ends of reads (see option -q in trim_galore manual) [25]
	-a adapter : (trim_galore) sequence of adapter for adapter trimming (default is universal Illumina adapter) ["AGATCGGAAGAG"]
	-m bismark_mismatch : (bismark) number of mismatches allowed in seed region [0]
	-l seedlen : (bismark) length of seed region [20]
	-L minIL : (STAR) minimum intron length (see STAR manual, --alignIntronMin option) [70]
	-I maxIL : (STAR) maximum intron length (see STAR manual, --alignIntronMax option) [5000]
	-N STAR_mismatch : (STAR) max number of mismatches, as a fraction of read length (for PE reads, use only length of single read since they're not treated as pairs) (see STAR -outFilterMismatchNoverReadLmax option) [0.05]
	-Q minmapQ : (STAR and bismark) min mapping quality (mapQ) in SAM alignment files for alignments to be retained; same value used for both STAR and bismark (see notes) [10]
	-x minreadsRNA : (scRNAseq analysis) min total RNA-seq read counts required to consider well valid [1000]
	-f mingenes : (scRNAseq analysis) min genes with at least one RNA-seq read count required to consider well valid [200]
	-X minbinsWGBS : (scBSseq analysis) min total fraction of bins (see -B) that must have at least read in order to consider well valid [0.1]
	-B binsize : (scBSseq analysis) calculate average methylation over bins this size (in bp) tiled genome-wide [1000]
	-T skipto : start analysis from this step (-G 1 will just start you at the beginning...) [1]
	-M memalloc : memory (in Gb) to allocate to java for running MarkDuplicates, make sure memalloc * threads is less than the total RAM available [5]
	-P picard : location of picard.jar from picard_tools (default assumes a variable named PICARD points to it) [$PICARD]
Flag options:
	-F : force allow reads from the two input files to not be same length (e.g if already trimmed; usually should not be used) [forcelen=false]
	-K : keep intermediate temporary files in $path_to_scratch (default these are deleted) [keeptemp=false]
	-r : allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!) [overwrite=false]
	-0 : checks that all required programs installed on PATH and all required helper scripts can be located, then exits without running
	-h : prints this version and usage information

Must be installed on your PATH:
	- bismark (by Babraham Bioinformatics, v.0.20.1)
	- bowtie2 (by Ben Langmead, v.2.3.4.3)
	- python (v3.9.6)
	- cutadapt (v.3.4 by Marcel Martin)
	- trim_galore (from Babraham Institute, tested on v0.6.7)
	- samtools (by Heng Li, tested on v1.11)
	- bedGraphToBigWig (part of UCSC tools, tested on v4)
	- htseq-count (part of HTSeq python package, tested on v0.13.5)
	- bedtools (v.2.30.0 by quinlanlab.org and others)
Must be in path_to_scripts - set path to this folder using -s option:
	- deduplicate_bismark.pl (by Babraham Bioinformatics, downloaded Dec. 18, 2018)
	- bismark_summary_to_bed.py (by Colette Picard, v.1.0, 10/26/2014)
	- steve_filter.py (v.1.0 by Colette Picard, 02/28/2020)
Other installation:
	- Markduplicates from picard_tools v.2.25.0 : set path to picard_tools jar using -P (default assumes it's saved to a variable named PICARD)
	
------------------------------------------------------------------------------------
EOF

[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# ----------------------
# Get user-specified arguments
# ----------------------

# Initiate environment
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 
workdir=$( pwd )											# working directory

# Required arguments:
# ----------------------
forward=""							# name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok)
bismark_index=""							# bismark-indexed reference genome for alignment (see bismark_genome_preparation)
STAR_index=""							# STAR-indexed transcriptome for alignment
gtf=""							# annotation file in GTF format
outdir=""							# name of directory with all output files; will be created if nonexistant

# Additional options:
# ----------------------
reverse=""							# name of file containing reverse reads if using paired-end data (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok)
name="expt"							# short experiment name used to name output files
path_to_scripts="$scriptDir/helper_scripts"							# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
path_to_scratch=""							# path to folder to use for scratch/fast calculations; WARNING: FILES AT THIS LOCATION MAY BE OVERWRITTEN, do not point to important folders (will be in $outdir by default, but you can speed up this script by pointing to a disk w/ faster I/O)
threads=1								# indicate how many threads are available to this program (note: bismark uses 4 threads by default, recommend using at least 4)
idxfile="$scriptDir/snmCAT_well_barcodes.txt"							# (demux) path to file containing well barcodes in fasta format
trim=10								# trims t bases from 5' and 3' end of reads before any additional quality- or adapter-based trimming (for read w/ well barcode, additional 8bp is added to this); recommend setting to 10 given adaptase behavior (see trim_galore option -t)
minQ=25								# (trim_galore) min quality to use for trimming off bad ends of reads (see option -q in trim_galore manual)
adapter="AGATCGGAAGAG"							# (trim_galore) sequence of adapter for adapter trimming (default is universal Illumina adapter)
bismark_mismatch=0								# (bismark) number of mismatches allowed in seed region
seedlen=20								# (bismark) length of seed region
minIL=70								# (STAR) minimum intron length (see STAR manual, --alignIntronMin option)
maxIL=5000								# (STAR) maximum intron length (see STAR manual, --alignIntronMax option)
STAR_mismatch=0.05								# (STAR) max number of mismatches, as a fraction of read length (for PE reads, use only length of single read since they're not treated as pairs) (see STAR -outFilterMismatchNoverReadLmax option)
minmapQ=10								# (STAR and bismark) min mapping quality (mapQ) in SAM alignment files for alignments to be retained; same value used for both STAR and bismark (see notes)
minreadsRNA=1000								# (scRNAseq analysis) min total RNA-seq read counts required to consider well valid
mingenes=200								# (scRNAseq analysis) min genes with at least one RNA-seq read count required to consider well valid
minbinsWGBS=0.1								# (scBSseq analysis) min total fraction of bins (see -B) that must have at least read in order to consider well valid
binsize=1000								# (scBSseq analysis) calculate average methylation over bins this size (in bp) tiled genome-wide
skipto=1								# start analysis from this step (-G 1 will just start you at the beginning...)
memalloc=5								# memory (in Gb) to allocate to java for running MarkDuplicates, make sure memalloc * threads is less than the total RAM available
picard="$PICARD"							# location of picard.jar from picard_tools (default assumes a variable named PICARD points to it)

# Flag options:
# ----------------------
forcelen=false							# force allow reads from the two input files to not be same length (e.g if already trimmed; usually should not be used)
keeptemp=false							# keep intermediate temporary files in $path_to_scratch (default these are deleted)
overwrite=false							# allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!)

checkdep=false

# ----------------------
while getopts "1:b:t:g:o:2:n:s:S:p:i:c:q:a:m:l:L:I:N:Q:x:f:X:B:T:M:P:FKr0h" opt; do
	case $opt in
		1)	# name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok)
			forward="$OPTARG"
			;;
		b)	# bismark-indexed reference genome for alignment (see bismark_genome_preparation)
			bismark_index="$OPTARG"
			;;
		t)	# STAR-indexed transcriptome for alignment
			STAR_index="$OPTARG"
			;;
		g)	# annotation file in GTF format
			gtf="$OPTARG"
			;;
		o)		# name of directory with all output files; will be created if nonexistant
			outdir="$OPTARG"
			;;
		2)	# name of file containing reverse reads if using paired-end data (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok)
			reverse="$OPTARG"
			;;
		n)	# short experiment name used to name output files
			name="$OPTARG"
			;;
		s)	# path to folder containing all required scripts (note - $scriptDir in default is the location of this script)
			path_to_scripts="$OPTARG"
			;;
		S)	# path to folder to use for scratch/fast calculations; WARNING: FILES AT THIS LOCATION MAY BE OVERWRITTEN, do not point to important folders (will be in $outdir by default, but you can speed up this script by pointing to a disk w/ faster I/O)
			path_to_scratch="$OPTARG"
			;;
		p)	# indicate how many threads are available to this program (note: bismark uses 4 threads by default, recommend using at least 4)
			threads="$OPTARG"
			;;
		i)	# (demux) path to file containing well barcodes in fasta format
			idxfile="$OPTARG"
			;;
		c)	# trims t bases from 5' and 3' end of reads before any additional quality- or adapter-based trimming (for read w/ well barcode, additional 8bp is added to this); recommend setting to 10 given adaptase behavior (see trim_galore option -t)
			trim="$OPTARG"
			;;
		q)	# (trim_galore) min quality to use for trimming off bad ends of reads (see option -q in trim_galore manual)
			minQ="$OPTARG"
			;;
		a)	# (trim_galore) sequence of adapter for adapter trimming (default is universal Illumina adapter)
			adapter="$OPTARG"
			;;
		m)	# (bismark) number of mismatches allowed in seed region
			bismark_mismatch="$OPTARG"
			;;
		l)	# (bismark) length of seed region
			seedlen="$OPTARG"
			;;
		L)	# (STAR) minimum intron length (see STAR manual, --alignIntronMin option)
			minIL="$OPTARG"
			;;
		I)	# (STAR) maximum intron length (see STAR manual, --alignIntronMax option)
			maxIL="$OPTARG"
			;;
		N)	# (STAR) max number of mismatches, as a fraction of read length (for PE reads, use only length of single read since they're not treated as pairs) (see STAR -outFilterMismatchNoverReadLmax option)
			STAR_mismatch="$OPTARG"
			;;
		Q)	# (STAR and bismark) min mapping quality (mapQ) in SAM alignment files for alignments to be retained; same value used for both STAR and bismark (see notes)
			minmapQ="$OPTARG"
			;;
		x)	# (scRNAseq analysis) min total RNA-seq read counts required to consider well valid
			minreadsRNA="$OPTARG"
			;;
		f)	# (scRNAseq analysis) min genes with at least one RNA-seq read count required to consider well valid
			mingenes="$OPTARG"
			;;
		X)	# (scBSseq analysis) min total fraction of bins (see -B) that must have at least read in order to consider well valid
			minbinsWGBS="$OPTARG"
			;;
		B)	# (scBSseq analysis) calculate average methylation over bins this size (in bp) tiled genome-wide
			binsize="$OPTARG"
			;;
		T)	# start analysis from this step (-G 1 will just start you at the beginning...)
			skipto="$OPTARG"
			;;
		M)	# memory (in Gb) to allocate to java for running MarkDuplicates, make sure memalloc * threads is less than the total RAM available
			memalloc="$OPTARG"
			;;
		P)	# location of picard.jar from picard_tools (default assumes a variable named PICARD points to it)
			picard="$OPTARG"
			;;
		F)	# force allow reads from the two input files to not be same length (e.g if already trimmed; usually should not be used)
			forcelen=true
			;;
		K)	# keep intermediate temporary files in $path_to_scratch (default these are deleted)
			keeptemp=true
			;;
		r)	# allow overwrite of existing files in provided outdir - default FALSE (WARNING: all files currently in outdir will be deleted!)
			overwrite=true
			;;
		0)	# check dependencies ok then exit
			checkdep=true
			;;
		h)	# print usage and version information to stdout and exit
			echo "$usage"
			exit 0
			;;
		\?)
			echo "Invalid option: -$OPTARG" >&2
			exit 1
			;;
		:)
			echo "Option -$OPTARG requires an argument." >&2
			exit 1
			;;
	esac
done

# Check that all programs required on PATH are installed
# ----------------------
command -v bismark >/dev/null 2>&1 || { echo "Error: bismark is required on PATH but was not found"; exit 1; }
command -v bowtie2 >/dev/null 2>&1 || { echo "Error: bowtie2 is required on PATH but was not found"; exit 1; }
command -v python >/dev/null 2>&1 || { echo "Error: python is required on PATH but was not found"; exit 1; }
command -v cutadapt >/dev/null 2>&1 || { echo "Error: cutadapt is required on PATH but was not found"; exit 1; }
command -v trim_galore >/dev/null 2>&1 || { echo "Error: trim_galore is required on PATH but was not found"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Error: samtools is required on PATH but was not found"; exit 1; }
command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "Error: bedGraphToBigWig is required on PATH but was not found"; exit 1; }
command -v htseq-count >/dev/null 2>&1 || { echo "Error: htseq-count is required on PATH but was not found"; exit 1; }
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools is required on PATH but was not found"; exit 1; }

# Check picard_tools is working
# ----------------------
res=$( java -jar "$picard" MarkDuplicates --version 2>&1 )
[[ "$res" =~ "Version"* ]] || { echo "Could not run MarkDuplicates from the picard-tools suite, check .jar lives in $picard (if not, change with -P); otherwise it may be incompatible with your java install OR you haven't allocated enough memory to java (see -m param, and check you have enough RAM)"; exit 1; }

# Check that required files are in path_to_scripts (set path to this location with option -s)
# ----------------------
[ ! -f "$path_to_scripts/deduplicate_bismark.pl" ] && { echo "Error: could not find required file deduplicate_bismark.pl in provided folder (${path_to_scripts})"; exit 1; }
[ ! -f "$path_to_scripts/bismark_summary_to_bed.py" ] && { echo "Error: could not find required file bismark_summary_to_bed.py in provided folder (${path_to_scripts})"; exit 1; }
[ ! -f "$path_to_scripts/steve_filter.py" ] && { echo "Error: could not find required file steve_filter.py in provided folder (${path_to_scripts})"; exit 1; }

# Done checking all requirements. Stop here if -0 flagged.
# ----------------------
"$checkdep" && exit 0

# Check all required inputs are provided
# ----------------------
[ -z "$forward" ] && { echo "Error: -1 forward is a required argument (name of file containing single-end reads or forward paired-end reads (FASTQ format, must have extension *.txt, *.fastq, or *.fq; if gzipped, then .txt.gz, .fastq.gz and .fq.gz all ok))"; exit 1; }
[ -z "$bismark_index" ] && { echo "Error: -b bismark_index is a required argument (bismark-indexed reference genome for alignment (see bismark_genome_preparation))"; exit 1; }
[ -z "$STAR_index" ] && { echo "Error: -t STAR_index is a required argument (STAR-indexed transcriptome for alignment)"; exit 1; }
[ -z "$gtf" ] && { echo "Error: -g gtf is a required argument (annotation file in GTF format)"; exit 1; }
[ -z "$outdir" ] && { echo "Error: -o outdir is a required argument (name of directory with all output files; will be created if nonexistant)"; exit 1; }

# Check that all required input files exist and can be accessed
# ----------------------
[ -f "$forward" ] || { echo "Could not open forward reads file -1 $forward"; exit 1; }
[[ ! -z "$reverse" && ! -f "$reverse" ]] && { echo "Could not open reverse reads file -2 $reverse"; exit 1; }
[ -f "$idxfile" ] || { echo "Could not open index file -i $idxfile"; exit 1; }
[ -f "$gtf" ] || { echo "Could not open GTF file -g $gtf"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.2.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.2.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.3.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.3.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.4.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.4.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.rev.1.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.rev.1.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.rev.2.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/CT_conversion/BS_CT.rev.2.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.2.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.2.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.3.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.3.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.4.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.4.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.rev.1.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.rev.1.bt2"; exit 1; }
[ -f "$bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.rev.2.bt2" ] || { echo "Could not open bismark_index file $bismark_index/Bisulfite_Genome/GA_conversion/BS_GA.rev.2.bt2"; exit 1; }
[ -f "$STAR_index/chrLength.txt" ] || { echo "Could not open STAR STAR_index file $STAR_index/chrLength.txt"; exit 1; }
[ -f "$STAR_index/chrName.txt" ] || { echo "Could not open STAR STAR_index file $STAR_index/chrName.txt"; exit 1; }
[ -f "$STAR_index/chrNameLength.txt" ] || { echo "Could not open STAR STAR_index file $STAR_index/chrNameLength.txt"; exit 1; }

# Check all other parameter values valid
# ----------------------
[[ "$threads" =~ ^[0-9]+$ ]] || { echo "Value of -p (num threads available) must be an integer; value provided was $threads"; exit 1; }
[[ "$trim" =~ ^[0-9]+$ ]] || { echo "Value of -c (num bases to trim from 3', 5' ends of reads) must be an integer; value provided was $trim"; exit 1; }
[[ "$minQ" =~ ^[0-9]+$ ]] || { echo "Value of -q (min base quality used for trimming) must be an integer; value provided was $minQ"; exit 1; }
[[ "$adapter" =~ ^[ATGC]+$ ]] || { echo "Value of -a (Illumina adapter to trim off reads) must contain only letters A, T, G and C; value provided was $adapter"; exit 1; }
[[ "$bismark_mismatch" =~ ^[0-9]+$ ]] || { echo "Value of -m (mismatches allowed for bismark) must be an integer; value provided was $bismark_mismatch"; exit 1; }
[[ "$seedlen" =~ ^[0-9]+$ ]] || { echo "Value of -l (length of seed region for bismark) must be an integer; value provided was $seedlen"; exit 1; }
[[ "$minIL" =~ ^[0-9]+$ ]] || { echo "Value of -L (min intron length for STAR) must be an integer; value provided was $minIL"; exit 1; }
[[ "$maxIL" =~ ^[0-9]+$ ]] || { echo "Value of -L (max intron length for STAR) must be an integer; value provided was $maxIL"; exit 1; }
[[ "$STAR_mismatch" =~ ^0.[0-9]+$ ]] || { echo "Value of -N (mismatches allowed for STAR) must be a fraction between [0,1) (example: 0.2, 0.567, etc.); value provided was $STAR_mismatch"; exit 1; }
[[ "$minmapQ" =~ ^[0-9]+$ ]] || { echo "Value of -Q (min mapping quality) must be an integer; value provided was $minmapQ"; exit 1; }
[[ "$skipto" =~ ^[0-9]+$ ]] || { echo "Value of -G (skip to this step) must be an integer; value provided was $skipto"; exit 1; }
[[ "$memalloc" =~ ^[0-9]+$ ]] || { echo "Value of -M (skip to this step) must be an integer; value provided was $memalloc"; exit 1; }
[[ "$minreadsRNA" =~ ^[0-9]+$ ]] || { echo "Value of -x (min total RNA-seq counts) must be an integer; value provided was $minreadsRNA"; exit 1; }
[[ "$minbinsWGBS" =~ ^0.[0-9]+$ ]] || { echo "Value of -X (min total WGBS reads) must be a value between [0,1); value provided was $minbinsWGBS"; exit 1; }
[[ "$binsize" =~ ^[0-9]+$ ]] || { echo "Value of -B (bin size for tiling genome) must be an integer; value provided was $binsize"; exit 1; }
[[ "$mingenes" =~ ^[0-9]+$ ]] || { echo "Value of -f (min RNA-seq detected genes) must be an integer; value provided was $mingenes"; exit 1; }

# If -G 1 and output folder doesn't exist, make it, else overwrite if user used -r, else error
# If -G >1 and output folder doesn't exist, error
# ----------------------
if [ "$skipto" -gt 1 ]; then
	[ ! -d "$outdir" ] && { echo "Could not find output directory ${outdir}; when using -G to skip earlier parts of this script, you must provide same -o outdir and -n name as original run."; exit 1; }
else
	if [ ! -d "$outdir" ]; then
		mkdir -p "$outdir"
	else
		if [ "$(ls -A $outdir)" ]; then
			if [ "$overwrite" = "true" ]; then
				echo "Overwriting previous contents of output dir $outdir"
				rm -rf "$outdir"/*	
			else
				echo "Error: provided output directory is not empty. To allow overwrite, use -r flag. WARNING: all existing files in -outdir- will be deleted. Seriously."
				exit 1
			fi
		fi
	fi
fi

# Make $SCRATCH folder, note that nothing is deleted but files may be overwritten since I don't check if it contains anything
if [ -z "$path_to_scratch" ]; then
	path_to_scratch="$outdir/scratch"
fi
if [ "$skipto" -gt 1 ]; then
	[ ! -d "$path_to_scratch" ] && { echo "Could not find scratch directory ${path_to_scratch}; when using -G to skip earlier parts of this script, you must provide same -o outdir, -S scratch and -n name as original run."; exit 1; }
else
	if [ ! -d "$path_to_scratch" ]; then
		mkdir -p "$path_to_scratch"
	else
		if [ "$(ls -A $path_to_scratch)" ]; then
			if [ "$overwrite" = "true" ]; then
				echo "Overwriting previous contents of scratch dir $path_to_scratch"
				rm -rf "$path_to_scratch"/*	
			else
				echo "Error: provided scratch directory $path_to_scratch is not empty. To allow overwrite, use -r flag. WARNING: all existing files in -path_to_scratch- will be deleted. Seriously."
				exit 1
			fi
		fi
	fi
fi

# Additional setup calculations
# ----------------------
# Each parallel instance of bismark takes ~4 threads, how many should we launch?
bismark_instances=$(( $threads / 4 ))
[[ "$bismark_instances" -eq 0 ]] && bismark_instances=1

# ----------------------
# Helper functions for this script
# ----------------------
err_msg ()
# prints an error message both to stdout and to the log file, then exits; usage: err_msg msg logfile
{
	echo "Error: $1"
	echo "Exited due to error: $1" >> $2
	exit 1	
}

displaytime () {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  [[ $D > 0 ]] && printf '%d days ' $D
  [[ $H > 0 ]] && printf '%d hours ' $H
  [[ $M > 0 ]] && printf '%d minutes ' $M
  [[ $D > 0 || $H > 0 || $M > 0 ]] && printf 'and '
  printf '%d seconds\n' $S
}

compress_sam_multi ()
# Small function to sort, compress and index (SAM->BAM) a list of SAM files (comma-separated, no spaces)
# sorts and compresses a SAM file to BAM and then indexes it if index=true (default false); 
# then removes the original SAM file. Each output is saved to the original filename (but .bam instead of .sam).
# Scriptlog is the log for this script.
# Usage: compress_sam_multi infiles
{
	index=false
	infile="$1"
	[ "$#" -eq 2 ] && index="$2"

	filelist=$( echo "$1" | tr ',' ' ' )

	# grab filename minus extension
	for infile in $filelist; do
		ff="${infile%.*}"
		if [[ -f "$infile" && -s "$infile" ]]; then
			samtools sort -@ "$threads" "$infile" -o "${ff}.bam" -T "${ff}" > "${ff}_log.txt" 2>&1
			[ $? != 0 ] && err_msg "samtools sort failed for input file $infile" "$log"
			if [ "$index" = "true" ]; then
				samtools index "${ff}.bam"
				[ $? != 0 ] && err_msg "samtools index failed for input file ${ff}.bam" "$log"
			fi
			rm "$infile" "${ff}_log.txt"
		else
			err_msg "tried to sort and index empty or nonexistant SAM file $infile" "$log"
		fi
	done
}

get_phred_qual ()
# Guess the quality encoding of the reads in a fastq file, by comparing against characters unique to phred+33 or phred+64
# Usage: get_phred_qual infile.fq
{
	infile="$1"
	phred33="\ ! \" # $ % & ' ( ) * + , - . / 0 1 2 3 4 5 6 7 8 9 : ; < = > ?"
	phred64="K L M N O P Q R S T U V W X Y Z [ \\ ] ^ _ \` a b c d e f g h i j k l m n o p q r s t u v w x y z { | } ~"

	# compare first 1000 reads against these chars
	tot33=0; tot64=0
	for char in $phred33; do
	res=$( awk 'NR % 4 == 0' $infile | head -1000 | awk -F"$char" '{tt+=NF-1} END {print tt}' )
	tot33=$(( $tot33 + $res ))
	done
	for char in $phred64; do
	res=$( awk 'NR % 4 == 0' $infile | head -1000 | awk -F"$char" '{tt+=NF-1} END {print tt}' )
	tot64=$(( $tot64 + $res ))
	done

	if [ "$tot33" -gt 0 ] && [ "$tot64" -gt 0 ]; then
		[ "$tot33" -gt "$tot64" ] && touse="phred33" || touse="phred64"
		echo "Warning: proper quality encoding could not be definitively determined. This shouldn't happen, so check your files. Using most likely encoding: $touse"
	elif [ "$tot33" -eq 0 ] && [ "$tot64" -eq 0 ]; then
		echo "Warning: input file $infile is not in proper fastq format (every 4th line should contain quality encodings, but this doesn't seem to be the case here)"
		touse="fail"
	else
		[ "$tot33" -gt 0 ] && touse="phred33" || touse="phred64"
	fi
	echo $touse
}

demux_PE ()
# Adds well-level info to read names in the provided pair of fastq files, combines output into single file
# and performs hard trimming of 5' and 3' ends based on $trim and desired final read length. Assumes well
# barcode is first 8bp of R1.
# Outputs three fastq files: a single file of successfully demuxed reads with R1 and R2 combined (_R1 and _R2 appended
# to read names), and an R1 and R2 file of reads that couldn't be assigned to a well barcode.
# Usage: demux_PE forward.fastq reverse.fastq barcodes.txt outprefix trim R1len R2len logfile
{
	if file "$1" | grep -q "gzip compressed"; then gz1="true"; else gz1="false"; fi
	if file "$2" | grep -q "gzip compressed"; then gz2="true"; else gz2="false"; fi
	awk -F$'\t' -v outp="$4" -v tt="$5" -v R1len="$6" -v R2len="$7" 'BEGIN {OFS=FS; totreads=0} 
	NR==FNR {
		idxmap[$2]=$1; readcounts[$2]=0; next
	} {
		if ($1 in idxmap) {
			printf("%s\n%s\n%s\n%s\n","@"idxmap[$1]"_R1_"substr($2,2),substr($3,tt+8+1,R1len),$4,substr($5,tt+8+1,R1len)) > outp"_barcoded.fastq"
			printf("%s\n%s\n%s\n%s\n","@"idxmap[$1]"_R2_"substr($6,2),substr($7,tt+1,R2len),$8,substr($9,tt+1,R2len)) > outp"_barcoded.fastq"
			readcounts[$1]+=1; totreads+=1
		} else {
			printf("%s\n%s\n%s\n%s\n",$2,$3,$4,$5) > outp"_demux_fail_R1.fastq"
			printf("%s\n%s\n%s\n%s\n",$6,$7,$8,$9) > outp"_demux_fail_R2.fastq"
			unassigned+=1; totreads+=1
		}
	} END {
		for (i in idxmap) {
			print idxmap[i],i,readcounts[i],(readcounts[i]/totreads)*100"%"
		}
		print "unassigned","N/A",unassigned,(unassigned/totreads)*100"%"
		print "total","N/A",totreads,"100%"
	}' "$3" <( paste <( sed 'N;N;N; s/\n/\t/g' <( [ "$gz1" = "true" ] && gunzip -c "$1" || cat "$1" ) ) <( sed 'N;N;N; s/\n/\t/g' <( [ "$gz2" = "true" ] && gunzip -c "$2" || cat "$2" ) ) | awk -F$'\t' '{OFS=FS} {idx=substr($2,1,8); if (idx !~ /N/) {print idx,$0}}' ) | sort -k1,1 > "$8"
}

demux_SE ()
# Same as demux_PE() but for single-end input
# Usage: demux_PE forward.fastq barcodes.txt outprefix trim R1len logfile
{
	if file "$1" | grep -q "gzip compressed"; then gz1="true"; else gz1="false"; fi
	awk -F$'\t' -v outp="$3" -v tt="$4" -v R1len="$5" 'BEGIN {OFS=FS; totreads=0} 
	NR==FNR {
		idxmap[$2]=$1; readcounts[$2]=0; next
	} {
		if ($1 in idxmap) {
			printf("%s\n%s\n%s\n%s\n","@"idxmap[$1]"_"substr($2,2),substr($3,tt+8+1,R1len),$4,substr($5,tt+8+1,R1len)) > outp"_barcoded.fastq"
			readcounts[$1]+=1; totreads+=1
		} else {
			printf("%s\n%s\n%s\n%s\n",$2,$3,$4,$5) > outp"_demux_fail.fastq"
			unassigned+=1; totreads+=1
		}
	} END {
		for (i in idxmap) {
			print idxmap[i],i,readcounts[i],(readcounts[i]/totreads)*100"%"
		}
		print "unassigned","N/A",unassigned,(unassigned/totreads)*100"%"
		print "total","N/A",totreads,"100%"
	}' "$2" <( sed 'N;N;N; s/\n/\t/g' <( [ "$gz1" = "true" ] && gunzip -c "$1" || cat "$1" ) | awk -F$'\t' '{OFS=FS} {idx=substr($2,1,8); if (idx !~ /N/) {print idx,$0}}' ) | sort -k1,1 > "$6"
}

dedup_RNA_multi ()
# Small function to run picard tools MarkDuplicates on multiple input files sequentially
# Input filelist should be comma-separated list (no spaces)
# Usage: dedup_RNA_multi filelist outprefix
# Example: dedup_RNA_multi file1.bam,file2.bam,file3.bam outprefix
# Output for each file is saved to filename_dedup.bam, log files are saved to outprefix_xyz
{
	filelist=$( echo "$1" | tr ',' ' ' ); outprefix="$2"
	for ff in $filelist; do
#		echo "Processing $ff"
		bb="${ff%.*}"; bname=$( basename $bb )
				
		java -Xmx${memalloc}g -jar $picard MarkDuplicates --I "$ff" --O "${bb}_dedup.bam" -METRICS_FILE "${outprefix}_${bname}_metrics.txt" -REMOVE_DUPLICATES true > "${outprefix}_${bname}_log.txt" 2>&1
		[ $? != 0 ] && err_msg "picard_tools MarkDuplicates failed for input file ${ff}, see log ${outprefix}_${bname}_log.txt" "$log"
		
		samtools index "${bb}_dedup.bam"
	done
	return 0
}

dedup_WGBS_multi ()
# Small function to run Bismark deduplicate_bismark.pl script on multiple input files sequentially
# Input filelist should be comma-separated list (no spaces)
# Usage: dedup_WGBS_multi filelist odir ldir scriptlog
# Example: dedup_WGBS_multi file1.bam,file2.bam,file3.bam odir ldir
# Output for each file is saved to filename_dedup.bam, log files are saved to odir_log_xyz
# odir = directory for output SAM files, ldir = directory for output log/report files
{
	filelist=$( echo "$1" | tr ',' ' ' ); odir="$2"; ldir="$3"
	for ff in $filelist; do
#		echo "Processing $ff"
		bb="${ff%.*}"; bname=$( basename $bb )
			
		$path_to_scripts/deduplicate_bismark.pl -s --output_dir "$odir" "$ff" > "$odir/dedup_log_${bname}.txt" 2>&1
		[ $? != 0 ] && err_msg "deduplicate_bismark.pl failed for input file ${ff}, see log $odir/dedup_log_${bname}.txt" "$log"
		rm "$odir/dedup_log_${bname}.txt"
		mv "$odir/${bname}.deduplication_report.txt" "$ldir/dedup_bismark_${bname}_summary.txt"
		mv "$odir/${bname}.deduplicated.sam" "$odir/${bname}_dedup.sam" 
	done
}

htseq_count_multi ()
# Small function to run htseq-count on multiple BAM files at once
# Usage: htseq_count_multi filelist odir
# filelist = comma-separated list of input BAM files (no spaces); odir = directory for output files
{
	filelist=$( echo "$1" | tr ',' ' ' ); odir="$2"
	for ff in $filelist; do
#		echo "htseq-count of $ff"
		well=$( echo "$ff" | sed -E "s/.+${name}_(.+)_RNAseq_dedup.bam/\1/" )	
		htseq-count -r pos --nonunique none -m intersection-nonempty -f bam -s no "$ff" "$gtf" > "$odir/${name}_${well}_counts.txt"	2> "$odir/${name}_${well}_log.txt"
		[ $? != 0 ] && err_msg "htseq-count failed for input file ${ff}, see $odir/${name}_${well}_log.txt" "$log"
		rm "$odir/${name}_${well}_log.txt"
	done
}

me_extractor_multi ()
# Small function to run bismark_methylation_extractor on multiple SAM files at once
# and convert to BED format, then get average across genome-wide tiles
# Saves output files to a separate folder for each input file, all in $odir
# Usage: bismark_methylation_extractor filelist odir_me odir_tiles
# filelist = comma-separated list of input SAM files (no spaces)
# odir_me = output directory for per-position methylation data
# odir_tiles = output directory for genome-wide tiling avg. methylation data
# uid = unique identifier for this batch (to avoid collisions w/ other instances of me_extractor_multi)
{
	filelist=$( echo "$1" | tr ',' ' ' ); odir="$2"; odir_avg="$3"; uid="$4"
	totbins=$( cat "$path_to_scratch/methylation_tiling/genome_windows.bed" | wc -l )

	# set up aggregate output files
	for context in CpG CHG CHH; do
		cat "$path_to_scratch/methylation_tiling/genome_windows_wheader.txt" > "$odir_avg/tmp${uid}_${context}_avg_all.txt"
		cat "$path_to_scratch/methylation_tiling/genome_windows_wheader.txt" > "$odir_avg/tmp${uid}_${context}_avg_QCpass.txt"
	done	
	printf "well\tme_CpG\tcov_CpG\tme_CHG\tcov_CHG\tme_CHH\tcov_CHH\tstatus\n" > "$odir_avg/tmp${uid}_summary.txt"
	
	# extract methylation from each sample
	wellno=1
	for ff in $filelist; do
#		echo "methylation_extractor of $ff"
		well=$( echo "$ff" | sed -E "s/.+${name}_(.+)_WGBS_dedup.sam/\1/" )	
		mkdir "$odir/${name}_${well}"; mkdir "$odir_avg/${name}_${well}"
		bismark_methylation_extractor -o "$odir/${name}_${well}" -s --no_header --report "$ff" 1> /dev/null 2> "$odir/${name}_${well}/log_methylation_extractor.txt"
		[ $? != 0 ] && err_msg "bismark_methylation_extractor failed for input file ${ff}, see $odir/${name}_${well}/log_methylation_extractor.txt" "$log"
		rm "$odir/${name}_${well}/log_methylation_extractor.txt"
		
		stat="fail"
		printf "$well" >> "$odir_avg/tmp${uid}_summary.txt"

		# convert results to bed format
		for context in "CpG" "CHG" "CHH"; do
			filelist=$( ls "$odir/${name}_${well}/${context}"* 2>/dev/null )
			fileno=$( echo $filelist | wc -w )
			if [ "$fileno" -gt 0 ]; then
				$path_to_scripts/bismark_summary_to_bed.py $filelist "$odir/${name}_${well}/${name}_${well}_${context}" > "$odir/${name}_${well}/log_bismark_summary_to_bed.txt"
				[ $? != 0 ] && err_msg "bismark_summary_to_bed.py failed for input file ${ff}, see $odir/${name}_${well}/log_bismark_summary_to_bed.txt" "$log"
				rm "$odir/${name}_${well}/log_bismark_summary_to_bed.txt" $filelist
				
				# sort resulting BED file
				sort -k1,1 -k2n,2 "$odir/${name}_${well}/${name}_${well}_${context}.bed" > "$odir/${name}_${well}/${name}_${well}_${context}_tmp.bed"
				mv "$odir/${name}_${well}/${name}_${well}_${context}_tmp.bed" "$odir/${name}_${well}/${name}_${well}_${context}.bed"
				
				# run bedtools map to get averages across tiles
				bedtools map -a "$path_to_scratch/methylation_tiling/genome_windows.bed" -b "$odir/${name}_${well}/${name}_${well}_${context}.bed" -c 4 -o sum | awk -F$'\t' '$4 != "."' > "$odir_avg/${name}_${well}/${name}_${well}_${context}_sum_unme.bed"
				bedtools map -a "$path_to_scratch/methylation_tiling/genome_windows.bed" -b "$odir/${name}_${well}/${name}_${well}_${context}.bed" -c 5 -o sum | awk -F$'\t' '$4 != "."' > "$odir_avg/${name}_${well}/${name}_${well}_${context}_sum_me.bed"
				
				# combine together
				bedtools intersect -a "$odir_avg/${name}_${well}/${name}_${well}_${context}_sum_unme.bed" -b "$odir_avg/${name}_${well}/${name}_${well}_${context}_sum_me.bed" -wa -wb | awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$4,$8,$8/($4+$8)}' > "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins.bed"
				rm "$odir_avg/${name}_${well}/${name}_${well}_${context}_sum_unme.bed" "$odir_avg/${name}_${well}/${name}_${well}_${context}_sum_me.bed"
				# make version w/ only bins with >= 5 total reads, drop if empty though
				awk '$4+$5 >= 5' "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins.bed" > "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed"
				[ ! -s "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed" ] && rm "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed"

				# summarize and add to aggregate files
				if [ -f "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed" ]; then
					nt=$( cat "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed" | wc -l ); ft=$(echo "scale=6; $nt / $totbins" | bc)
					avgme=$( cut -f6 "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed" | awk '{ss+=$1; tt+=1} END {print ss/tt * 100}' )
					[ $( echo "$ft > $minbinsWGBS" | bc -l ) -eq 1 ] && stat="pass"
					printf "\t%0.3f%%\t%0.2f%%" "$(echo "scale=6; $ft * 100" | bc)" "$avgme" >> "$odir_avg/tmp${uid}_summary.txt"
				else
					printf "\t0\tN/A" >> "$odir_avg/tmp${uid}_summary.txt"
				fi
			else
				printf "\t0\tN/A" >> "$odir_avg/tmp${uid}_summary.txt"	
			fi		

			# add to full dataset
			curpos=$(( $wellno + 2 )); newpos=$(( $curpos + 6 ))
			echo "$well" > "$odir_avg/tmp${uid}_tmp_tiles.bed"
			if [ -f "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed" ]; then
				# merge back to full list of tiles
				bedtools intersect -a "$path_to_scratch/methylation_tiling/genome_windows.bed" -b "$odir_avg/${name}_${well}/${name}_${well}_${context}_over_${binsize}bp_bins_min5.bed" -loj | cut -f9 | awk -F$'\t' '{OFS=FS} $NF=="." {$NF=""}1' >> "$odir_avg/tmp${uid}_tmp_tiles.bed"
			else
				yes "0" | head -n $totbins >> "$odir_avg/tmp${uid}_tmp_tiles.bed"
			fi
			paste "$odir_avg/tmp${uid}_${context}_avg_all.txt" "$odir_avg/tmp${uid}_tmp_tiles.bed" > "$odir_avg/tmp${uid}_tmp.txt"
			mv "$odir_avg/tmp${uid}_tmp.txt" "$odir_avg/tmp${uid}_${context}_avg_all.txt"
			if [ "$stat" = "pass" ]; then
				paste "$odir_avg/tmp${uid}_${context}_avg_QCpass.txt" "$odir_avg/tmp${uid}_tmp_tiles.bed" > "$odir_avg/tmp${uid}_tmp.txt"
				mv "$odir_avg/tmp${uid}_tmp.txt" "$odir_avg/tmp${uid}_${context}_avg_QCpass.txt"
			fi
		done
		wellno=$(( wellno + 1 ))
		printf "\t$stat\n" >> "$odir_avg/tmp${uid}_summary.txt"				
	done
	rm -f "$odir_avg/tmp${uid}_tmp_tiles.bed" 
}


# ----------------------
# Main code
# ----------------------
log="$outdir/${name}_log.txt" 	# create log file
mkdir -p "$outdir/logs"			# folder to store log files
mkdir -p "$outdir/summaries"			# folder to store summary files
time_start=$(date)				# time run was started
time_ss=$(date +%s)				# time run was started (in seconds)

# Output user-derived options to stdout and to log file
# ----------------------

# If script was run before, check command is the same
if [ "$skipto" -gt 1 ]; then
[ -f "$log" ] || { echo "Could not locate log file (expected ${log}); make sure when running with -G that the rest of the command is unchanged."; exit 1; }
# TODO
echo "" | tee -a "$log"
echo "Running snmCAT_seq_map v1.0 by Colette Picard (03/15/2022):" | tee -a "$log"
echo "Run start on: $time_start" | tee -a "$log"
echo "Starting analysis at step: $skipto" | tee -a "$log"

# Else if this is first time running script, write summary of all parameters:
else
echo "Running snmCAT_seq_map v1.0 by Colette Picard (03/15/2022):" | tee "$log"
echo "Run start on: $time_start" | tee -a "$log"
fi
echo "-------------------------" | tee -a "$log"
echo "Working directory: $( pwd )" | tee -a "$log"
[ -z "$reverse" ] && echo "Reads file: $forward" | tee -a "$log"
[ -z "$reverse" ] || echo "Forward reads file: $forward" | tee -a "$log"
[ -z "$reverse" ] || echo "Reverse reads file: $reverse" | tee -a "$log"
echo "Output directory: $outdir" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Bismark index: $bismark_index" | tee -a "$log"
echo "STAR index: $STAR_index" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Additional settings:" | tee -a "$log"
echo "Stubname: $name" | tee -a "$log"
echo "Helper scripts are in: $path_to_scripts" | tee -a "$log"
echo "Scratch space for temp I/O: $path_to_scratch" | tee -a "$log"
echo "Number of threads available: $threads" | tee -a "$log"
echo "Plate well index file (fasta): $idxfile" | tee -a "$log"
echo "Read trimming settings (trim_galore):" | tee -a "$log"
echo "   # bases to trim at 5' and 3' ends of each read (excluding the 8bp well index): $trim" | tee -a "$log"
echo "   Quality cutoff for trimming at 3' end: $minQ" | tee -a "$log"
echo "   Adapter sequence to trim: $adapter" | tee -a "$log"
echo "Bisulfite-seq alignment settings (bismark):" | tee -a "$log"
echo "   Number of instances of bismark to launch (~4 threads each): $bismark_instances" | tee -a "$log"
echo "   Length of seed region: $seedlen" | tee -a "$log"
echo "   Number of mismatches allowed in seed: $bismark_mismatch" | tee -a "$log"
echo "RNA-seq alignment settings (STAR):" | tee -a "$log"
echo "   Minimum intron length: $minIL" | tee -a "$log"
echo "   Maximum intron length: $maxIL" | tee -a "$log"
echo "   Max number of mismatches as a fraction of total read length: $STAR_mismatch" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Full command used:" | tee -a "$log"
echo "$0 $@" | tee -a "$log"
echo "-------------------------" | tee -a "$log"


# ----------------------
# Step 0: check all inputs ok
# ----------------------
printf "\nStep 0: Checking that all inputs are ok\n" | tee -a "$log"

# Get all chromosomes in provided file, check bismark and STAR indexes have same chromosomes
# ----------------------
allchrs=$( bowtie2-inspect "$bismark_index/Bisulfite_Genome/CT_conversion/BS_CT" -n | cut -f1 -d' ' | sed 's/_CT_converted//' | tr '\n' ' ' | tr '\r' ' ' | tr '\r\n' ' ' )
allchrsSTAR=$( cat "$STAR_index/chrName.txt" )
allchrsarray=( $allchrs )
allchrsarraySTAR=( $allchrsSTAR )
res=$( echo ${allchrsarray[@]} ${allchrsarraySTAR[@]} | tr ' ' '\n' | sort | uniq -c | sed 's/ \+/\t/g' | cut -f2 | sort -u )			# == 2 if all chrs match
[ "$res" -ne 2 ] && err_msg "chromosomes in STAR index and Bismark index don't match" "$log"
numchrs="${#allchrsarray[@]}"
echo "$numchrs chromosomes detected in genome indexes: $allchrs" | tee -a "$log"

# Check that input FASTQ file(s)s have correct extensions
# ----------------------
gzipped=false
f_basename=$( basename $forward ); ext_f="${f_basename#*.}"
f_dirname=$(dirname $forward)
[ ! -z "$reverse" ] && { r_basename=$( basename $reverse ); ext_r="${r_basename#*.}"; }
[ ! -z "$reverse" ] && r_dirname=$(dirname $reverse)
[[ ! -z "$reverse" && "$ext_f" != "$ext_r" ]] && err_msg "Input files $forward and $reverse must have the same file extension (because this saves the script author some time); note this error may also occur if you have dots (.) in your filenames. Don't put dots in your filenames except to indicate file extensions (general life advice haha)." "$log"

if [ "$ext_f" != "txt" ] && [ "$ext_f" != "fq" ] && [ "$ext_f" != "fastq" ]; then
	if [ "$ext_f" = "txt.gz" ] || [ "$ext_f" = "fq.gz" ] || [ "$ext_f" = "fastq.gz" ]; then
		gzipped=true
	else
		err_msg "Input file(s) -1, -2 have unrecognized extension $ext_f , permitted extensions are .txt, .fq, and .fastq (must be FASTQ format) or .txt.gz, .fq.gz, .fastq.gz if gzipped; note this error may also occur if you have dots (.) in your filenames. Don't put dots in your filenames except to indicate file extensions (general life advice haha)." "$log"
	fi
fi

# If extension indicates gzipped file, check that input files are actually gzipped; if not, check -not- gzipped
# ----------------------
if "$gzipped"; then
	gz1=false; gz2=false
	if file "$forward" | grep -q "gzip compressed"; then gz1=true; echo "Input file $forward is compressed" ; else err_msg "Input file $forward has extension ending in .gz but does not appear to be a gzipped file" "$log"; fi
	[ ! -z "$reverse" ] && { if file "$reverse" | grep -q "gzip compressed"; then gz2=true; echo "Input file $reverse is compressed" ; else err_msg "Input file $reverse has extension ending in .gz but does not appear to be a gzipped file" "$log"; fi; }
else
	if file "$forward" | grep -q "gzip compressed"; then err_msg "Input file $forward does not have extension ending in .gz but does appear to be a gzipped file" "$log" ; else echo "File $forward is not compressed"; fi
	[ ! -z "$reverse" ] && { if file "$reverse" | grep -q "gzip compressed"; then err_msg "Input file $reverse does not have extension ending in .gz but does appear to be a gzipped file" "$log" ; else echo "File $reverse is not compressed"; fi; }
fi

# Get the quality encoding of the reads
# ----------------------
if "$gzipped"; then
	gunzip -c "$forward" | head -4000 > "$outdir/tmp_f.fq"
	[ ! -z "$reverse" ] && { gunzip -c $reverse | head -4000 > "$outdir/tmp_r.fq"; }
else
	head -4000 "$forward" > "$outdir/tmp_f.fq"
	[ ! -z "$reverse" ] && { head -4000 $reverse > "$outdir/tmp_r.fq"; }
fi

if [ ! -z "$reverse" ]; then
	qualR1=$( get_phred_qual "$outdir/tmp_f.fq" )
	qualR2=$( get_phred_qual "$outdir/tmp_r.fq" )
	[ "$qualR1" = "fail" ] && { echo "Error: could not determine quality encoding of ${forward}, is it in proper fastq format?"; exit 1; }
	[ "$qualR2" = "fail" ] && { echo "Error: could not determine quality encoding of ${reverse}, is it in proper fastq format?"; exit 1; }
	[ "$qualR1" = "$qualR2" ] || { echo "Error: forward and reverse reads file do not have same quality encoding (forward is ${qualR1}, reverse is ${qualR2})"; exit 1; }
	res=$qualR1
else
	res=$( get_phred_qual "$outdir/tmp_f.fq" )
	[ "$qualR1" = "fail" ] && { echo "Error: could not determine quality encoding of ${forward}, is it in proper fastq format?"; exit 1; }
fi

rm -f "$outdir/tmp_f.fq" "$outdir/tmp_r.fq"
[ "$res" = "phred33" ] && phred33=true || phred33=false
[ "$res" = "phred33" ] && { echo "Quality encoding was auto-detected as PHRED+33" | tee -a "$log";} || { echo "Quality encoding was auto-detected as PHRED+64" | tee -a "$log"; }

# Check expected index sequences can be detected in R1 (this is a quick method so doesn't fully check, just checks first 10,000 reads)
# ----------------------
if [ "$skipto" -eq 1 ]; then
	# reformat barcodes file
	sed '$!N;s/\n/\t/' "$idxfile" | sed 's/[>^]//g' > "$path_to_scratch/barcodes_reformatted.txt"
	# grab barcodes of first 10000 reads to test if they're in the file
	if "$gzipped"; then	
		gunzip -c "$forward" | head -40000 | awk 'NR % 4 == 2' | cut -c-8 | awk '$0 !~ /N/' > "$path_to_scratch/barcodes_test.txt"
	else
		cat "$forward" | head -40000 | awk 'NR % 4 == 2' | cut -c-8 | awk '$0 !~ /N/' > "$path_to_scratch/barcodes_test.txt"
	fi
	# how many barcodes are in the index file?
	res=$( awk 'NR==FNR {idx[$2]=$1;next} idx[$1]' $path_to_scratch/barcodes_reformatted.txt $path_to_scratch/barcodes_test.txt | wc -l )
	[ "$res" -le 1000 ] && err_msg "in small initial check, fewer than 10% of reads could be assigned to a barcode in index file ${idxfile}, please check that file is correct" "$log"
	[ "$res" -le 5000 ] && echo "Warning: in small initial check, fewer than 50% of reads could be assigned to a barcode in index file ${idxfile}, please check that file is correct"
	rm "$path_to_scratch/barcodes_test.txt"
fi

welllist=$( cut -f1 "$path_to_scratch/barcodes_reformatted.txt" )		# use this a bit later


# Get initial number of reads in input file(s); if paired-end then check that # reads same in both files
# ----------------------
if [ "$skipto" -eq 1 ]; then
	if "$gzipped"; then
		f_tot=$( gunzip -c "$forward" | wc -l )								# of lines in forward file
		[ ! -z "$reverse" ] && r_tot=$( gunzip -c "$reverse" | wc -l )		# of lines in reverse file
		readf=$( gunzip -c "$forward" | head -n2 | tail -1 ); readlenf=$( echo ${#readf} )
		[ ! -z "$reverse" ] && { readr=$( gunzip -c "$reverse" | head -n2 | tail -1 ); readlenr=$( echo ${#readr} ); }
	else
		f_tot=$( wc -l "$forward" | awk '{print $1}' )								# of lines in forward file
		[ ! -z "$reverse" ] && r_tot=$( wc -l "$reverse" | awk '{print $1}' )		# of lines in reverse file
		readf=$( head -n2 $forward | tail -1 ); readlenf=$( echo ${#readf} )
		[ ! -z "$reverse" ] && { readr=$( head -n2 $reverse | tail -1 ); readlenr=$( echo ${#readr} ); }
	fi

	[ $(( $f_tot % 4 )) -ne 0 ] && err_msg "number of lines in -1 fastq file ($f_tot) is not a multiple of four (expected for fastq file)." "$log"
	[[ ! -z "$reverse" && $(( $r_tot % 4 )) -ne 0 ]] && err_msg "number of lines in -2 fastq file ($r_tot) is not a multiple of four (expected for fastq file)." "$log"

	# Get # of reads by dividing by 4, check these are same
	f_tot=$(( $f_tot/4 ))
	[ ! -z "$reverse" ] && r_tot=$(( $r_tot/4 ))
	[[ ! -z "$reverse" && "$f_tot" -ne "$r_tot" ]] && err_msg "-1 and -2 files do not contain same number of reads.\n # Forward reads = $f_tot,  # Reverse reads = $r_tot" "$log"

	# also get read lengths (assume all reads same length, so use length of read in first entry)
	if [ "$forcelen" = "false" ]; then
		[[ ! -z "$reverse" && "$readlenf" -ne "$readlenr" ]] && err_msg "Reads in -1 and -2 are not the same length (based on first read only); override with -F" "$log"
	fi
	printf "\nStarting number of reads: %s" "$f_tot" | tee -a "$log"	
	[ ! -z "$reverse" ] && printf "\nReads are ${readlenf}x${readlenr}bp\n" | tee -a "$log"	
	[ ! -z "$reverse" ] || printf "\nReads are ${readlenf}bp\n" | tee -a "$log"	
	echo "-------------------------" | tee -a "$log"
fi


# ----------------------
# Step 1: demux reads according to well
# ----------------------
if [ "$skipto" -le 1 ]; then
	ts=$(date +%s)
	printf "\nStep 1: demultiplexing reads according to plate well barcode...\n" | tee -a "$log"
	mkdir -p "$path_to_scratch/demux_fasta"; mkdir -p "$outdir/fastqc_prefiltering"

	# -----------RUN FASTQC ON ALL READS-----------
	echo " - Running initial QC check with fastqc..." | tee -a "$log"
	fastqc "$forward" -t "$threads" -o "$outdir/fastqc_prefiltering" 2> /dev/null
	[ ! -z "$reverse" ] && fastqc "$reverse" -t "$threads" -o "$outdir/fastqc_prefiltering" 2> /dev/null

	# -----------DEMUX READS-----------
	oop="$path_to_scratch/demux_fasta/${name}_barcoded.fastq"
	oof="$path_to_scratch/demux_fasta/${name}_no_barcode"
	
	# note - to keep things simple, all this does is:
	# (1) add the well barcode to the read name (e.g. @F11_A00817:359:H7NHHDSX3:4:1101:1588:1016 1:N:0:CGGAAGATAA+AATATGCCAG_R1 -> well is F11)
	# (2) combine R1 and R2 into one file if PE data provided
	# (3) filter out reads w/o match to barcodes
	# (4) perform hard trimming of 5' and 3' ends of reads (see -t param, for R1 5' end adds 8bp to trim off well barcode too
	# Actual demux is done post-mapping since there's not much point doing it so early (and dealing with
	# so many small files is a bit of a pain).
	
	# For hard trimming, awk substr(str, start, len) function takes start and len arguments, so get
	# what these values should be for R1 (and R2 if relevant)
	R1len=$( echo "$readlenf - 2*${trim} - 8" | bc )
	R2len=$( echo "$readlenf - 2*${trim}" | bc )
	
	# 'Demultiplex' files, taking advantage of number of threads by splitting input into $threads subprocesses
	numl=$(( $(( $(( $f_tot / $threads )) + 1 )) * 4 ))
	echo " - Demultiplexing reads in $threads parallel instances of $numl reads each..." | tee -a "$log"
	[ "$gz1" = "true" ] && gunzip -c "$forward" | split -d -a 3 -l $numl - "$path_to_scratch/demux_fasta/${name}_R1_tmp" || split -a 3 -d -l $numl "$forward" "$path_to_scratch/demux_fasta/${name}_R1_tmp"
	[ ! -z "$reverse" ] && { [ "$gz2" = "true" ] && gunzip -c "$reverse" | split -d -a 3 -l $numl - "$path_to_scratch/demux_fasta/${name}_R2_tmp" || split -a 3 -d -l $numl "$reverse" "$path_to_scratch/demux_fasta/${name}_R2_tmp"; }

	pids=(); cmds=()
	for ffR1 in "$path_to_scratch/demux_fasta/${name}_R1_"*; do
		fileno=$( echo "$ffR1" | sed -E "s/.+R1_tmp(.+)/\1/" )
		if [ ! -z "$reverse" ]; then
			ffR2=$( echo "$ffR1" | sed "s/_R1_tmp/_R2_tmp/" )
#			echo "Demuxing ${ffR1} and ${ffR2}..."
			cmd="demux_PE $ffR1 $ffR2 $path_to_scratch/barcodes_reformatted.txt $path_to_scratch/demux_fasta/${name}_${fileno} $trim $R1len $R2len $outdir/logs/demux_summary_initial_${fileno}.txt"
			eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
		else
#			echo "Demuxing ${ffR1}..."
			cmd="demux_SE $ffR1 $path_to_scratch/barcodes_reformatted.txt $path_to_scratch/demux_fasta/${name}_${fileno} $trim $R1len $outdir/logs/demux_summary_initial_${fileno}.txt"
			eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
		fi	
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "demultiplexing failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# merge back all files and delete temp files
	echo " - All $threads complete, merging files..." | tee -a "$log"
	rm "$path_to_scratch/demux_fasta/${name}_R"*	
#	cat $path_to_scratch/demux_fasta/*_barcoded.fastq > "$path_to_scratch/demux_fasta/${name}_barcoded_all.fastq"
#	rm $path_to_scratch/demux_fasta/*_barcoded.fastq	
	if [ ! -z "$reverse" ]; then
		cat $path_to_scratch/demux_fasta/*_fail_R1.fastq | gzip > $path_to_scratch/demux_fasta/${name}_demux_failed_R1.fastq
		cat $path_to_scratch/demux_fasta/*_fail_R2.fastq | gzip > $path_to_scratch/demux_fasta/${name}_demux_failed_R2.fastq
		rm $path_to_scratch/demux_fasta/*_fail_*.fastq
	else
		cat $path_to_scratch/demux_fasta/*_fail.fastq | gzip > $path_to_scratch/demux_fasta/${name}_demux_failed.fastq
		rm $path_to_scratch/demux_fasta/*_fail.fastq
	fi	
	printf "well\tindex\treads_sequenced\tpct_of_pool\n" > "$outdir/summaries/demux_summary_initial.txt"
	awk -F$'\t' '{OFS=FS} {counts[$1"\t"$2]+=$3} END {for (aa in counts) {print aa,counts[aa],counts[aa]/counts["total\tN/A"]*100"%"}}' $outdir/logs/demux_summary_initial_*.txt | sort -k1,1 >> "$outdir/summaries/demux_summary_initial.txt"
	[ "$keeptemp" = "false" ] && rm $outdir/logs/demux_summary_initial_*.txt
		
	# Summarize demux results
	n_tot=$( awk '/total/ {print $3}' "$outdir/summaries/demux_summary_initial.txt" )
	n_wel=$( head -n -2 "$outdir/summaries/demux_summary_initial.txt" | wc -l ); well50=$(( n_wel / 2 )); well25=$(( n_wel / 4 )); well75=$(( n_wel - well25 ))
	n_unk=$( awk '/unassigned/ {print $3}' "$outdir/summaries/demux_summary_initial.txt" ); p_unk=$(echo "scale=3; $n_unk / $n_tot *100" | bc)
	n_kno=$(( $n_tot - $n_unk )); p_kno=$(echo "scale=3; $n_kno / $n_tot *100" | bc)
	n_25=$( cut -f3 "$outdir/summaries/demux_summary_initial.txt" | head -n -2 | sort -k1n,1 | head -${well25} | tail -1 ); p_25=$(echo "scale=3; $n_25 / $n_tot *100" | bc)
	n_med=$( cut -f3 "$outdir/summaries/demux_summary_initial.txt" | head -n -2 | sort -k1n,1 | head -${well50} | tail -1 ); p_med=$(echo "scale=3; $n_med / $n_tot *100" | bc)
	n_75=$( cut -f3 "$outdir/summaries/demux_summary_initial.txt" | head -n -2 | sort -k1n,1 | head -${well75} | tail -1 ); p_75=$(echo "scale=3; $n_75 / $n_tot *100" | bc)
	n_avg=$( cut -f3 "$outdir/summaries/demux_summary_initial.txt" | head -n -2 | awk '{s+=$1} END {print s/NR}' ); p_avg=$(echo "scale=3; $n_avg / $n_tot *100" | bc)

	[ -z "$reverse" ] && readstr="reads" || readstr="read pairs"
	printf " - Summary of demultiplexing results (note read pairs count as 1 read):\n" | tee -a "$log"
	printf "   - Number of $readstr processed: $n_tot \n" | tee -a "$log"
	printf "   - Successfully assigned to wells: $n_kno (%0.1f%%)\n" "$p_kno" | tee -a "$log" 
	printf "   - Unassigned: $n_unk (%0.1f%%)\n" "$p_unk" | tee -a "$log" 
	printf "   - Stats for number of $readstr assigned to each barcode:\n" | tee -a "$log" 
	printf "     - Average $readstr per well: $n_avg (%0.1f%%)\n" "$p_avg" | tee -a "$log" 
	printf "     - 25th percentile $readstr per well: $n_25 (%0.1f%%)\n" "$p_25" | tee -a "$log" 
	printf "     - Median $readstr per well: $n_med (%0.1f%%)\n" "$p_med" | tee -a "$log" 
	printf "     - 75th percentile $readstr per well: $n_75 (%0.1f%%)\n" "$p_75" | tee -a "$log" 
	echo " - Full demux results summary saved to: $outdir/summaries/demux_summary_initial.txt" | tee -a "$log"
	
	[ -z "$reverse" ] && echo "** Note: each read of a pair is counted as two separate reads for the rest of this script, since pairs are considered separately **"
	
	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

else
	# if skipping this step, make sure required output files exist
	printf "\nSkipping step 1. Checking all step 1 required outputs exist...\n" | tee -a "$log"
	if [ "$skipto" -eq 2 ]; then
		[ -f "$path_to_scratch/demux_fasta/${name}_000_barcoded.fastq" ] || err_msg "option -T was used to skip step 1, but could not find file generated in step 1 required to run step 2 (${path_to_scratch}/demux_fasta/${name}_000_barcoded.fastq)" "$log"
	fi
fi


# ----------------------
# Step 2: filter and trim reads using trim_galore
# ----------------------
if [ "$skipto" -le 2 ]; then
	ts=$(date +%s)
	printf "\nStep 2: filtering and trimming reads using trim_galore\n" | tee -a "$log"
	mkdir -p "$path_to_scratch/qc_filtered"; mkdir -p "$outdir/fastqc_postfiltering"

	# at this point, reads have been pseudo-demuxed (well info is first few chars of read name)
	# and if PE, first and second read files have been combined

	[ "$phred33" = "true" ] && qualstr="phred33" || qualstr="phred64"

	# -----------RUN TRIM_GALORE-----------

	# run this as $threads subprocesses
	echo " - Trimming and filtering low-quality reads in $threads parallel instances..." | tee -a "$log"
	pids=(); cmds=()
	for ff in "$path_to_scratch/demux_fasta/${name}"*"barcoded.fastq"; do
		fileno=$( echo "$ff" | sed -E "s/.+_([0-9]+)_barcoded.+/\1/" )
		cmd="trim_galore --fastqc --fastqc_args \"--outdir $outdir/fastqc_postfiltering\" --${qualstr} -a $adapter -q $minQ --stringency 3 -o $path_to_scratch/qc_filtered/${name}_${fileno} $ff > $outdir/logs/log_trim_galore_${fileno}.txt 2>&1"
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done

	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "trim_galore failed, command that failed: ${cmds[i]}" "$log"
	done

	# summarize results
	totQC=0; removedQC=0
	for ll in "$outdir/logs/log_trim_galore_"*; do
		totQC=$( awk -v tp=$"totQC" '/sequences processed in total/ {print $1+tp}' "$ll" )
		removedQC=$( awk -v tf=$"removedQC" '/Sequences removed/ {print $14+tf}' "$ll" )
	done
	remQC=$(( $totQC - $removedQC ))
	[ "$remQC" = "0" ] && err_msg "No reads passed quality filtering step, please check your filtering parameters (in particular, was PHRED+ encoding detected incorrectly? Please report bugs!)" "$log"
	echo " - Read trimming and filtering complete, outputting summary:" | tee -a "$log"
	echo "   - Number of reads analyzed: $totQC" | tee -a "$log"
	printf "   - Number reads removed (read < 20bp): $removedQC (%0.1f%%)\n" "$(echo "scale=3; $removedQC / $totQC *100" | bc)" | tee -a "$log"
	printf "   - Number reads remaining: $remQC (%0.1f%%)\n" "$(echo "scale=3; $remQC / $totQC *100" | bc)" | tee -a "$log"
	rm "$outdir/logs/log_trim_galore_"*

	echo " - Merging outputs and deleting temp files..." | tee -a "$log"
	cat "$path_to_scratch/qc_filtered/${name}_"*"/"*"trimmed"* > "$path_to_scratch/qc_filtered/${name}_filt.fastq"
	rm -rf "$path_to_scratch/qc_filtered/${name}_"*/
	
	# if everything above successful, delete temp files...
	[ "$keeptemp" = "false" ] && rm -rf "$path_to_scratch/demux_fasta"

	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	
else
	# if skipping this step, make sure required output files exist
	printf "Skipping step 2. Checking all step 2 required outputs exist...\n" | tee -a "$log"
	if [ "$skipto" -eq 3 ]; then
		[ -f "$path_to_scratch/qc_filtered/${name}_filt.fastq" ] || err_msg "option -T was used to skip step 2, but could not find file generated in step 2 required for later steps ($path_to_scratch/qc_filtered/${name}_filt.fastq)" "$log"
	fi
fi


# ----------------------
# Step 3: align reads to transcriptome using STAR
# ----------------------
if [ "$skipto" -le 3 ]; then
	ts=$(date +%s)
	printf "\nStep 3: aligning reads to transcriptome using STAR\n" | tee -a "$log"
	
	mkdir -p "$path_to_scratch/STAR"
	
	# -----------RUN STAR-----------
	echo " - Performing RNA-seq alignment with STAR using $threads threads..." | tee -a "$log"
	STAR --alignEndsType EndToEnd --runThreadN "$threads" --genomeDir "$STAR_index" --readFilesIn "$path_to_scratch/qc_filtered/${name}_filt.fastq" --outFilterType BySJout --outFilterMultimapNmax 1 --winAnchorMultimapNmax 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax "$STAR_mismatch" --alignIntronMin "$minIL" --alignIntronMax "$maxIL" --outFileNamePrefix "$path_to_scratch/STAR/${name}" --outSAMtype BAM Unsorted --outFilterIntronMotifs RemoveNoncanonical --outReadsUnmapped Fastx > "$outdir/logs/STAR_log.txt"
	[ $? != 0 ] && err_msg "error when aligning reads using STAR, see $outdir/logs/STAR_log.txt" "$log"
	
	# rename files and delete temp files
	mv "$path_to_scratch/STAR/${name}Log.final.out" "$outdir/logs/STAR_result_summary.txt"
	rm "$path_to_scratch/STAR/${name}Log.out" "$path_to_scratch/STAR/${name}Log.progress.out"
	mv "$path_to_scratch/STAR/${name}SJ.out.tab" "$path_to_scratch/STAR/all_detected_junctions.txt"
	mv "$path_to_scratch/STAR/${name}Aligned.out.bam" "$path_to_scratch/STAR/STAR_alignments_raw.bam"
	mv "$path_to_scratch/STAR/${name}Unmapped.out.mate1" "$path_to_scratch/STAR/STAR_unmapped.fastq"
	rm -rf "$path_to_scratch/STAR/${name}_STARtmp"
	
	# filter out low mapping quality reads and convert to SAM ** note: w/o header! **
	samtools view -@ "$threads" -q "$minmapQ" "$path_to_scratch/STAR/STAR_alignments_raw.bam" > "$path_to_scratch/STAR/STAR_alignments.sam"
	[ $? != 0 ] && err_msg "samtools view failed, command: samtools view -@ $threads -hq $minmapQ $path_to_scratch/bismark/bismark_alignments_raw.bam > $path_to_scratch/bismark/bismark_alignments.sam" "$log"
	postfilt=$( cat "$path_to_scratch/STAR/STAR_alignments.sam" | wc -l )
	
	# done, summarize results
	echo " - Done aligning to STAR, brief alignment summary:" | tee -a "$log"
	totin=$( awk '/Number of input reads/ {print $6}' "$outdir/logs/STAR_result_summary.txt" )
	printf "   - A total of $postfilt (%0.1f%%) reads out of $totin mapped uniquely to transcriptome with mapQ > ${minmapQ}\n" "$(echo "scale=3; $postfilt / $totin *100" | bc)" | tee -a "$log"
	
	# convert to SAM and demux
	echo " - Demultiplexing STAR alignment files..." | tee -a "$log"
	# add SAM header to all files first
	for well in $welllist; do
		[ "${#well}" -eq 2 ] && touse="${well}_" || touse="$well"
		samtools view -H "$path_to_scratch/STAR/STAR_alignments_raw.bam" > "$path_to_scratch/STAR/${name}_${touse}_RNAseq.sam"
	done
	
	# demux SAM file
	printf "well\treads_obtained\n" > "$outdir/summaries/demux_summary_RNAseq.txt"
	cat "$path_to_scratch/STAR/STAR_alignments.sam" | awk -F$'\t' -v oo="$path_to_scratch/STAR/${name}" '{OFS=FS} {counts[substr($1,1,3)]+=1; print $0 >> oo"_"substr($1,1,3)"_RNAseq.sam"} END {for (a in counts) {print a,counts[a]}}' | sed 's/_//' | sort -k1,1 > "$outdir/summaries/demux_summary_RNAseq.txt"

	# minor renaming to deal w/ those double underscores... also delete empty files
	for ff in "$path_to_scratch/STAR/${name}_"*"_RNAseq.sam"; do
		rr=$( tail "$ff" | awk '$0 !~ /@/' | wc -l )
		if [ "$rr" -eq 0 ]; then
			well=$( echo "$ff" | sed -E "s/.*${name}_(.*)_.*/\1/" | sed 's/_//' )
			echo "Warning: no RNA-seq reads detected for well ${well}"
			rm "$ff"
		else	
			orgname=$( basename $ff )
			newname=$( echo $orgname | sed 's/__/_/' )
			ppath=$( dirname $ff )
			[ "$newname" != "$orgname" ] && mv "$ff" "$ppath/$newname"
		fi
	done
	
	# sort and compress all files in batches of (# of files/wells / # of threads)
	echo " - Sorting and compressing all files..." | tee -a "$log"
	
	numwells=$( ls -1 "$path_to_scratch/STAR/"*"_RNAseq.sam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/STAR/"*"_RNAseq.sam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/STARlist"
	
	pids=(); cmds=()
	for ff in "$path_to_scratch/STARlist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="compress_sam_multi $filelist" 
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "SAM file sorting failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# clean up temp files
	rm "$path_to_scratch/STARlist"*
	[ "$keeptemp" = "false" ] && rm -rf "$path_to_scratch/qc_filtered"
	
	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	
else
	printf "Skipping step 3. Checking all step 3 required outputs exist..." | tee -a "$log"
	rnaseqfiles=$( ls -1 "$path_to_scratch/STAR/${name}_"*"_RNAseq.bam" 2> /dev/null | wc -l )
	[ "$rnaseqfiles" -eq 0 ] && err_msg "option -T was used to skip step 3, but could not find any RNA-seq alignment files generated in step 3 (expected in $path_to_scratch/STAR)" "$log"
	echo "detected RNA-seq alignment files for $rnaseqfiles wells"
fi


# ----------------------
# Step 4: align reads to bisulfite-treated genome using bismark
# ----------------------
if [ "$skipto" -le 4 ]; then
	ts=$(date +%s)
	printf "\nStep 4: aligning remaining non-RNA-seq reads to bisulfite-treated genome using bismark\n" | tee -a "$log"
	
	mkdir -p "$path_to_scratch/bismark"; mkdir -p "$path_to_scratch/bismark/tmp"

	[ "$phred33" = "true" ] && qualstr="--phred33-quals" || qualstr="--phred64-quals"
	echo " - Performing WGBS alignment with bismark, launching $bismark_instances separate instances (~4 threads each)..." | tee -a "$log"

	# -----------RUN BISMARK-----------
	bismark "$qualstr" -N "$bismark_mismatch" -L "$seedlen" --non_directional -o "$path_to_scratch/bismark" --parallel "$bismark_instances" --temp_dir "$path_to_scratch/bismark/tmp" "$bismark_index" "$path_to_scratch/STAR/STAR_unmapped.fastq" > "$outdir/logs/bismark_log.txt" 2>&1
	[ $? != 0 ] && err_msg "error when aligning reads using bismark, see $outdir/logs/bismark_log.txt" "$log"
	
	# rename files and delete tmp files
	rm -rf "$path_to_scratch/bismark/tmp"
	mv "$path_to_scratch/bismark/"*"_report.txt" "$outdir/logs/bismark_result_summary.txt"
	mv "$path_to_scratch/bismark/"*".bam" "$path_to_scratch/bismark/bismark_alignments_raw.bam"

	# filter out low mapping quality reads and convert to SAM ** note: w/o header! **
	samtools view -@ "$threads" -q "$minmapQ" "$path_to_scratch/bismark/bismark_alignments_raw.bam" > "$path_to_scratch/bismark/bismark_alignments.sam"
	[ $? != 0 ] && err_msg "samtools view failed, command: samtools view -@ $threads -hq $minmapQ $path_to_scratch/bismark/bismark_alignments_raw.bam > $path_to_scratch/bismark/bismark_alignments.sam" "$log"
	postfilt=$( cat "$path_to_scratch/bismark/bismark_alignments.sam" | wc -l )

	# done, summarize results
	echo " - Done aligning to bismark, brief alignment summary:" | tee -a "$log"
	totin=$( awk '/Sequences analysed in total/ {print $5}' "$outdir/logs/bismark_result_summary.txt" )
	printf "   - A total of $postfilt (%0.1f%%) reads out of $totin mapped uniquely to the bisulfite-converted genome with mapQ > ${minmapQ}\n" "$(echo "scale=3; $postfilt / $totin *100" | bc)" | tee -a "$log"

	# filter out nonconverted and any leftover RNA-seq reads using the Steve filter script
	echo " - Filtering out nonconverted and any leftover RNA-seq reads using Steve filter with $threads instances" | tee -a "$log"	
	# split input file into $threads chunks of numl lines
	numl=$(( $(( $postfilt / $threads )) + 1 ))
	split -d -l $numl "$path_to_scratch/bismark/bismark_alignments.sam" "$path_to_scratch/bismark/bismark_alignments_tmp"

	pids=(); cmds=()
	for ff in "$path_to_scratch/bismark/bismark_alignments_tmp"*; do
		cmd="$path_to_scripts/steve_filter.py $ff $path_to_scratch/bismark/bismark_alignments_stevefilt_tmp_${ff: -2} > $path_to_scratch/bismark/stevefilt_${ff: -2}_log.txt"
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "steve_filter.py failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# successful, delete temp files, recombine files and summarize
	remreads=0
	for ff in "$path_to_scratch/bismark/bismark_alignments_tmp"*; do
		rem=$( awk '/alignments remain/ {print $2}' "$path_to_scratch/bismark/stevefilt_${ff: -2}_log.txt" )
		remreads=$(( remreads + rem ))
	done	
	cat "$path_to_scratch/bismark/bismark_alignments_stevefilt_tmp_"*"_pass.sam" > "$path_to_scratch/bismark/bismark_alignments_stevefilt_pass.sam"
	cat "$path_to_scratch/bismark/bismark_alignments_stevefilt_tmp_"*"_fail.sam" > "$path_to_scratch/bismark/bismark_alignments_stevefilt_fail.sam"
	rm "$path_to_scratch/bismark/bismark_alignments_tmp"* "$path_to_scratch/bismark/bismark_alignments_stevefilt_tmp_"*	
	printf "   - $remreads (%0.1f%%) reads out of $postfilt remain after Steve filter\n" "$(echo "scale=3; $remreads / $postfilt *100" | bc)" | tee -a "$log"
	
	# convert to SAM and demux
	echo " - Demultiplexing bismark alignment files..." | tee -a "$log"
	
	# add SAM header to all files first
	for well in $welllist; do
		[ "${#well}" -eq 2 ] && touse="${well}_" || touse="$well"
		samtools view -H "$path_to_scratch/bismark/bismark_alignments_raw.bam" > "$path_to_scratch/bismark/${name}_${touse}_WGBS.sam"
	done
	
	# demux SAM file
	printf "well\treads_obtained\n" > "$outdir/summaries/demux_summary_WGBS.txt"
	cat "$path_to_scratch/bismark/bismark_alignments_stevefilt_pass.sam" | awk -F$'\t' -v oo="$path_to_scratch/bismark/${name}" '{OFS=FS} {counts[substr($1,1,3)]+=1; print $0 >> oo"_"substr($1,1,3)"_WGBS.sam"} END {for (a in counts) {print a,counts[a]}}' | sed 's/_//' | sort -k1,1 >> "$outdir/summaries/demux_summary_WGBS.txt"

	# minor renaming to deal w/ those double underscores...
	for ff in "$path_to_scratch/bismark/${name}_"*"_WGBS.sam"; do
		rr=$( tail "$ff" | awk '$0 !~ /@/' | wc -l )
		if [ "$rr" -eq 0 ]; then
			well=$( echo "$ff" | sed -E "s/.*${name}_(.*)_.*/\1/" | sed 's/_//' )
			echo "Warning: no WGBS reads detected for well ${well}"
			rm "$ff"
		else	
			orgname=$( basename $ff )
			newname=$( echo $orgname | sed 's/__/_/' )
			ppath=$( dirname $ff )
			[ "$newname" != "$orgname" ] && mv "$ff" "$ppath/$newname"
		fi
	done
	
	# sort and compress all files in batches of (# of files/wells / # of threads)
	echo " - Sorting and compressing all files..." | tee -a "$log"
	
	numwells=$( ls -1 "$path_to_scratch/bismark/"*"_WGBS.sam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/bismark/"*"_WGBS.sam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/bismarklist"
	
	pids=(); cmds=()
	for ff in "$path_to_scratch/bismarklist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="compress_sam_multi $filelist" 
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "SAM file sorting failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# clean up temp files
	rm "$path_to_scratch/bismarklist"*

	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	
else
	printf "Skipping step 4. Checking all step 4 required outputs exist..." | tee -a "$log"
	bismarkfiles=$( ls -1 "$path_to_scratch/bismark/${name}_"*"_WGBS.bam" 2> /dev/null | wc -l )
	[ "$bismarkfiles" -eq 0 ] && err_msg "option -T was used to skip step 4, but could not find any WGBS alignment files generated in step 4 (expected in $path_to_scratch/bismark)" "$log"
	echo "detected WGBS alignment files for $bismarkfiles wells"
fi


# ----------------------
# Step 5: checking library complexity and removing PCR duplicates
# ----------------------
if [ "$skipto" -le 5 ]; then
	ts=$(date +%s)
	printf "\nStep 5: checking library complexity and removing PCR duplicates\n" | tee -a "$log"
	
	mkdir -p "$path_to_scratch/dedup_logs"
	
	numwells=$( ls -1 "$path_to_scratch/STAR/"*"_RNAseq.bam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/STAR/"*"_RNAseq.bam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/STARlist"
	echo " - Removing PCR duplicates from RNA-seq alignments in $threads batches of $wellsperthread wells..." | tee -a "$log"
	
	pids=(); cmds=()
	for ff in "$path_to_scratch/STARlist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="dedup_RNA_multi $filelist $path_to_scratch/dedup_logs/MarkDuplicates_log"
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "deduplication of RNA-seq data failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# successful, delete temp files and summarize results
	echo " - Done, summary of % duplicates rate across wells:" | tee -a "$log"
	rm -f "$path_to_scratch/STARlist"* "$path_to_scratch/dedup_logs/MarkDuplicates"*"log.txt"
	printf "well\ttotreads\tduplicates\tremaining\tpct_dup\n" > "$outdir/summaries/dedup_summary_RNAseq.txt"
	
	for ff in "$path_to_scratch/dedup_logs/MarkDuplicates"*; do
		wellno=$( echo "$ff" | sed -E "s/.+MarkDuplicates_log_${name}_(.+)_RNAseq_metrics.txt/\1/" )
		totreads=$( tail -3 "$ff" | head -1 | cut -f2 )
		dups=$( tail -3 "$ff" | head -1 | cut -f6 )
		rem=$(( totreads - dups ))
		pdup=$( echo "scale=5; $dups / $totreads *100" | bc )
		printf "${wellno}\t${totreads}\t${dups}\t${rem}\t%0.2f%%\n" "$pdup" >> "$outdir/summaries/dedup_summary_RNAseq.txt"
	done
			
	echo "   - average: "$( cut -f5 "$outdir/summaries/dedup_summary_RNAseq.txt" | tail -n+2 | awk '{s+=$1} END {printf("%0.2f%%",s/NR)}' ) | tee -a "$log"
	echo "   - 25th percentile: "$( cut -f5 "$outdir/summaries/dedup_summary_RNAseq.txt" | tail -n+2 | sort -k1n,1 | head -$(( numwells / 4 )) | tail -1 ) | tee -a "$log"
	echo "   - median: "$( cut -f5 "$outdir/summaries/dedup_summary_RNAseq.txt" | tail -n+2 | sort -k1n,1 | head -$(( numwells / 2 )) | tail -1 ) | tee -a "$log"
	echo "   - 75th percentile: "$( cut -f5 "$outdir/summaries/dedup_summary_RNAseq.txt" | tail -n+2 | sort -k1n,1 | head -$(( 3 * $(( numwells / 4 )) )) | tail -1 ) | tee -a "$log"

	# Repeat for WGBS data
	numwells=$( ls -1 "$path_to_scratch/bismark/"*"_WGBS.bam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/bismark/"*"_WGBS.bam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/bismarklist"
	echo " - Removing PCR duplicates from WGBS alignments in $threads batches of $wellsperthread wells..." | tee -a "$log"
	
	pids=(); cmds=()
	for ff in "$path_to_scratch/bismarklist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="dedup_WGBS_multi $filelist $path_to_scratch/bismark $path_to_scratch/dedup_logs"
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "deduplication of WGBS data failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# successful, summarize results	
	echo " - Done, summary of % duplicates rate across wells:" | tee -a "$log"
	rm -f "$path_to_scratch/bismarklist"*
	printf "well\ttotreads\tduplicates\tremaining\tpct_dup\n" > "$outdir/summaries/dedup_summary_WGBS.txt"

	for ff in "$path_to_scratch/dedup_logs/dedup_bismark"*; do
		wellno=$( echo "$ff" | sed -E "s/.+dedup_bismark_${name}_(.+)_WGBS_summary.txt/\1/" )
		totreads=$( awk '/Total number of alignments analysed/ {print $8}' "$ff" )
		dups=$( awk '/Total number duplicated alignments/ {print $6}' "$ff" )
		rem=$(( totreads - dups ))
		pdup=$( echo "scale=5; $dups / $totreads *100" | bc )
		printf "${wellno}\t${totreads}\t${dups}\t${rem}\t%0.2f%%\n" "$pdup" >> "$outdir/summaries/dedup_summary_WGBS.txt"
	done
	
	echo "   - average: "$( cut -f5 "$outdir/summaries/dedup_summary_WGBS.txt" | tail -n+2 | awk '{s+=$1} END {printf("%0.2f%%",s/NR)}' ) | tee -a "$log"
	echo "   - 25th percentile: "$( cut -f5 "$outdir/summaries/dedup_summary_WGBS.txt" | tail -n+2 | sort -k1n,1 | head -$(( numwells / 4 )) | tail -1 ) | tee -a "$log"
	echo "   - median: "$( cut -f5 "$outdir/summaries/dedup_summary_WGBS.txt" | tail -n+2 | sort -k1n,1 | head -$(( numwells / 2 )) | tail -1 ) | tee -a "$log"
	echo "   - 75th percentile: "$( cut -f5 "$outdir/summaries/dedup_summary_WGBS.txt" | tail -n+2 | sort -k1n,1 | head -$(( 3 * $(( numwells / 4 )) )) | tail -1 ) | tee -a "$log"

	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

else
	printf "Skipping step 5. Checking all step 5 required outputs exist...\n" | tee -a "$log"
	bismarkfiles=$( ls -1 "$path_to_scratch/bismark/${name}_"*"_WGBS_dedup.sam" 2> /dev/null | wc -l )
	rnaseqfiles=$( ls -1 "$path_to_scratch/STAR/${name}_"*"RNAseq_dedup.bam" 2> /dev/null | wc -l )
	[ "$bismarkfiles" -eq 0 ] && err_msg "option -T was used to skip step 5, but could not find any WGBS alignment post-PCR-dedup files generated in step 5 (expected in $path_to_scratch/bismark)" "$log"
	[ "$rnaseqfiles" -eq 0 ] && err_msg "option -T was used to skip step 5, but could not find any RNA-seq alignment post-PCR-dedup files generated in step 5 (expected in $path_to_scratch/STAR)" "$log"
	echo " - detected post-PCR-dedup RNA-seq alignment files for $rnaseqfiles wells"
	echo " - detected post-PCR-dedup WGBS alignment files for $bismarkfiles wells"
fi


# ----------------------
# Step 6: post-processing of RNA-seq data (read counting, etc.)
# ----------------------
if [ "$skipto" -le 6 ]; then
	ts=$(date +%s)
	printf "\nStep 6: getting count matrices from RNA-seq data\n" | tee -a "$log"

	# for each well, count reads over genes using htseq-count
	mkdir -p "$path_to_scratch/htseq_count"
	echo " - Getting read counts over features in $gtf using htseq-count" | tee -a "$log"
	numwells=$( ls -1 "$path_to_scratch/STAR/"*"_RNAseq_dedup.bam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/STAR/"*"_RNAseq_dedup.bam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/STARlist"
		
	pids=(); cmds=()
	for ff in "$path_to_scratch/STARlist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="htseq_count_multi $filelist $path_to_scratch/htseq_count"	
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "htseq-count failed, command that failed: ${cmds[i]}" "$log"
	done

	rm -f "$path_to_scratch/STARlist"*

	# combine all counts files into a big matrix file, one with all 384 wells regardless of status, one only with 'good' wells
	# include summary of # genes detected + average coverage of detected genes
	echo " - Done, combining into count matrices and outputting summary" | tee -a "$log"
	mkdir -p "$outdir/RNAseq_count_matrices"
	
	echo "gene_id" > "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt"
	cut -f1 $( ls "$path_to_scratch/htseq_count/"*"_counts.txt" | head -1 ) | awk '$0 !~ /^__/' >> "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt"
	cat "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt" > "$outdir/RNAseq_count_matrices/RNAseq_counts_QCpass.txt"
	
	printf "well\tgenes_detected\ttot_counts\tavg_per_gene\tno_feature\tstatus\n" > "$outdir/summaries/RNAseq_counts_summary.txt"
	passQC=0; numwells=$( echo $welllist | wc -w ); numgenes=$( tail -n+2 "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt" | wc -l )
	for ww in $welllist; do
		ff="$path_to_scratch/htseq_count/${name}_${ww}_counts.txt"
		if [ -f "$ff" ]; then
			stat="fail"
			# get number of genes detected, total cov. across those genes, and # reads assigned to no feature
			detected=$( awk -F$'\t' 'BEGIN {OFS=FS; ss=0} $1 !~ /^__/ && $2 > 0 {ss+=1} END {print ss}' "$ff" )
			sum_detected=$( awk -F$'\t' 'BEGIN {OFS=FS; ss=0} $1 !~ /^__/ && $2 > 0 {ss+=$2} END {print ss}' "$ff" )
			no_feature=$( grep "__no_feature" "$ff" | cut -f2 )
			[ "$detected" -gt 0 ] && avg_per_gene=$( echo "scale=8; $sum_detected / $detected" | bc ) || avg_per_gene="0"
			
			# add this well to *all* matrix
			echo "$ww" > "$path_to_scratch/htseq_count/tmp.txt"
			awk -F$'\t' '$1 !~ /^__/' | cut -f2 "$ff" >> "$path_to_scratch/htseq_count/tmp.txt"
			paste "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt" "$path_to_scratch/htseq_count/tmp.txt" > "$path_to_scratch/htseq_count/tmp2.txt"
			mv "$path_to_scratch/htseq_count/tmp2.txt" "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt"

			# does this well pass basic QC? if yes add to the passQC matrix as well
			if [ "$detected" -ge "$mingenes" ] && [ "$sum_detected" -ge "$minreadsRNA" ]; then
				paste "$outdir/RNAseq_count_matrices/RNAseq_counts_QCpass.txt" "$path_to_scratch/htseq_count/tmp.txt" > "$path_to_scratch/htseq_count/tmp2.txt"
				mv "$path_to_scratch/htseq_count/tmp2.txt" "$outdir/RNAseq_count_matrices/RNAseq_counts_QCpass.txt"
				passQC=$(( passQC + 1 )); stat="pass"
			fi
			
			printf "$ww\t$detected\t$sum_detected\t%0.1f\t$no_feature\t$stat\n" "$avg_per_gene" >> "$outdir/summaries/RNAseq_counts_summary.txt"
		else
			echo "$ww" > "$path_to_scratch/htseq_count/tmp.txt"
			yes "0" | head -n $numgenes >> "$path_to_scratch/htseq_count/tmp.txt"
			paste "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt" "$path_to_scratch/htseq_count/tmp.txt" > "$path_to_scratch/htseq_count/tmp2.txt"
			mv "$path_to_scratch/htseq_count/tmp2.txt" "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt"		
			printf "$ww\t0\t0\t0\t0\tfail\n" >> "$outdir/summaries/RNAseq_counts_summary.txt"
		fi
	done
	rm "$path_to_scratch/htseq_count/tmp.txt"
	avg_genes_pass=$( tail -n+2 "$outdir/summaries/RNAseq_counts_summary.txt" | awk 'BEGIN {tt=0; ss=0} $6=="pass" {ss+=$2; tt+=1} END {if (tt == 0) {print 0} else {print ss/tt}}' )
	avg_counts_pass=$( tail -n+2 "$outdir/summaries/RNAseq_counts_summary.txt" | awk 'BEGIN {tt=0; ss=0} $6=="pass" {ss+=$3; tt+=1} END {if (tt == 0) {print 0} else {print ss/tt}}' )	
	printf "   - Out of $numwells wells, $passQC (%0.1f%%) passed basic QC filters (at least $minreadsRNA total counts over $mingenes or more genes)\n" "$(echo "scale=3; $passQC / $numwells *100" | bc)" | tee -a "$log"
	printf "     - Avg. number of genes detected in wells that passed QC: %0.1f\n" "$avg_genes_pass" | tee -a "$log"
	printf "     - Avg. total counts obtained in wells that passed QC: %0.1f\n" "$avg_counts_pass" | tee -a "$log"

	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

else
	printf "Skipping step 6. Checking all step 6 required outputs exist..." | tee -a "$log"
	[ -f "$outdir/RNAseq_count_matrices/RNAseq_counts_all.txt" ] || err_msg "option -T was used to skip step 6, but could not find count matrices generated in step 6 (expected $outdir/RNAseq_count_matrices/RNAseq_counts_all.txt)" "$log"
fi


# ----------------------
# Step 7: post-processing of WGBS data (methylation extractor, getting average counts over genome tiles)
# ----------------------
if [ "$skipto" -le 7 ]; then
	ts=$(date +%s)
	printf "\nStep 7: extracting methylation info from WGBS data\n" | tee -a "$log"

	mkdir -p "$path_to_scratch/methylation_extractor"
	mkdir -p "$path_to_scratch/methylation_tiling"
	mkdir -p "$outdir/WGBS_tiling_matrices"
	# get chromosome sizes file from a BAM header
	samtools view -H "$path_to_scratch/bismark/bismark_alignments_raw.bam" | awk '$1 ~ /@SQ/ {print $2,$3}' | sed 's/[SL]N://g' | tr ' ' '\t' > "$outdir/chrom.sizes"
	bedtools makewindows -g "$outdir/chrom.sizes" -w "$binsize" > "$path_to_scratch/methylation_tiling/genome_windows.bed"
	printf "chr\tstart\tend\n" > "$path_to_scratch/methylation_tiling/genome_windows_wheader.txt"
	cat "$path_to_scratch/methylation_tiling/genome_windows.bed" >> "$path_to_scratch/methylation_tiling/genome_windows_wheader.txt"
	
	# process all samples in parallel depending on threads
	echo " - Extracting methylation data and getting avg. (weighted) methylation over ${binsize}bp bins tiled genome-wide" | tee -a "$log"
	numwells=$( ls -1 "$path_to_scratch/bismark/"*"_WGBS_dedup.sam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/bismark/"*"_WGBS_dedup.sam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/bismarklist"
		
	pids=(); cmds=(); uid=1
	for ff in "$path_to_scratch/bismarklist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="me_extractor_multi $filelist $path_to_scratch/methylation_extractor $path_to_scratch/methylation_tiling $uid"
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" ); uid=$(( uid + 1 ))
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "htseq-count failed, command that failed: ${cmds[i]}" "$log"
	done

	rm -f "$path_to_scratch/bismarklist"*
	
	# done, clean up, combine files and output summary
	echo " - Done, combining into tiling matrices and outputting summary" | tee -a "$log"
	for context in CpG CHG CHH; do
		# combine tiling matrix files
		for ff in "$path_to_scratch/methylation_tiling/tmp"*"_${context}_avg_"*".txt"; do
			bb="${ff%.*}"
			cut -f4- "$ff" > "${bb}_tmp.txt"
		done
		paste "$path_to_scratch/methylation_tiling/genome_windows_wheader.txt" "$path_to_scratch/methylation_tiling/tmp"*"_${context}_avg_all_tmp.txt" > "$outdir/WGBS_tiling_matrices/WGBS_${context}me_over_${binsize}bp_bins_all.txt"
		paste "$path_to_scratch/methylation_tiling/genome_windows_wheader.txt" "$path_to_scratch/methylation_tiling/tmp"*"_${context}_avg_QCpass_tmp.txt" > "$outdir/WGBS_tiling_matrices/WGBS_${context}me_over_${binsize}bp_bins_QCpass.txt"
		rm "$path_to_scratch/methylation_tiling/tmp"*"_${context}_avg_all_tmp.txt"
		rm "$path_to_scratch/methylation_tiling/tmp"*"_${context}_avg_QCpass_tmp.txt"

		# combine summary files
		cat "$path_to_scratch/methylation_tiling/tmp"*"_summary.txt" | awk 'NR!=1 && $0 !~ /cov.+status/' | sort -k1,1 > "$outdir/summaries/WGBS_methylation_summary.txt"
	done
		
	# output brief summary
	echo " - $numwells wells had aligning WGBS reads. Average stats for all wells:" | tee -a "$log"
	avgCGme=$( cut -f3 "$outdir/summaries/WGBS_methylation_summary.txt" | tail -n+2 | awk '{ss+=$1; tt+=1} END {print ss/tt}' )
	avgCGcov=$( cut -f2 "$outdir/summaries/WGBS_methylation_summary.txt" | tail -n+2 | awk '{ss+=$1; tt+=1} END {print ss/tt}' )
	avgCHGme=$( cut -f5 "$outdir/summaries/WGBS_methylation_summary.txt" | tail -n+2 | awk '{ss+=$1; tt+=1} END {print ss/tt}' )
	avgCHGcov=$( cut -f4 "$outdir/summaries/WGBS_methylation_summary.txt" | tail -n+2 | awk '{ss+=$1; tt+=1} END {print ss/tt}' )
	avgCHHme=$( cut -f7 "$outdir/summaries/WGBS_methylation_summary.txt" | tail -n+2 | awk '{ss+=$1; tt+=1} END {print ss/tt}' )
	avgCHHcov=$( cut -f6 "$outdir/summaries/WGBS_methylation_summary.txt" | tail -n+2 | awk '{ss+=$1; tt+=1} END {print ss/tt}' )		
	printf "   - CpG: had coverage over %0.3f%% of genome with avg. methylation %0.2f%%\n" "$avgCGcov" "$avgCGme" | tee -a "$log"
	printf "   - CHG: had coverage over %0.3f%% of genome with avg. methylation %0.2f%%\n" "$avgCHGcov" "$avgCHGme" | tee -a "$log"
	printf "   - CHH: had coverage over %0.3f%% of genome with avg. methylation %0.2f%%\n" "$avgCHHcov" "$avgCHHme" | tee -a "$log"
	
	passQC=$( cut -f8 "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $1=="pass" {tt+=1} END {print tt}' )
	printf " - Of these, $passQC wells passed QC filters (>= $minbinsWGBS frac of genome covered). " | tee -a "$log"
	if [ "$passQC" -gt 0 ]; then
		echo "Average stats for wells that passed QC:" | tee -a "$log"
		avgCGme=$( cat "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $8=="pass" {ss+=$3; tt+=1} END {if (tt==0) {print 0} else {print ss/tt}}' )
		avgCGcov=$( cat "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $8=="pass" {ss+=$2; tt+=1} END {if (tt==0) {print 0} else {print ss/tt}}' )
		avgCHGme=$( cat "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $8=="pass" {ss+=$5; tt+=1} END {if (tt==0) {print 0} else {print ss/tt}}' )
		avgCHGcov=$( cat "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $8=="pass" {ss+=$4; tt+=1} END {if (tt==0) {print 0} else {print ss/tt}}' )
		avgCHHme=$( cat "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $8=="pass" {ss+=$7; tt+=1} END {if (tt==0) {print 0} else {print ss/tt}}' )
		avgCHHcov=$( cat "$outdir/summaries/WGBS_methylation_summary.txt" | awk 'BEGIN{tt=0} $8=="pass" {ss+=$6; tt+=1} END {if (tt==0) {print 0} else {print ss/tt}}' )		
		printf "   - CpG: had coverage over %0.3f%% of genome with avg. methylation %0.2f%%\n" "$avgCGcov" "$avgCGme" | tee -a "$log"
		printf "   - CHG: had coverage over %0.3f%% of genome with avg. methylation %0.2f%%\n" "$avgCHGcov" "$avgCHGme" | tee -a "$log"
		printf "   - CHH: had coverage over %0.3f%% of genome with avg. methylation %0.2f%%\n" "$avgCHHcov" "$avgCHHme" | tee -a "$log"
	else
		echo "" | tee -a "$log"
	fi
				
	te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	
else
	printf "\nSkipping step 7. Checking all step 7 required outputs exist..." | tee -a "$log"
	for context in CpG CHG CHH; do
		[ -f "$outdir/WGBS_tiling_matrices/WGBS_${context}me_over_${binsize}bp_bins_all.txt" ] || err_msg "option -T was used to skip step 7, but could not find tiling methylation matrices generated in step 7 (expected $outdir/WGBS_tiling_matrices/WGBS_${context}me_over_${binsize}bp_bins_all.txt)" "$log"
	done
fi


# ----------------------
# Step 8: final summary and cleanup
# ----------------------
if [ "$skipto" -le 8 ]; then
	ts=$(date +%s)
	printf "\nStep 8: performing final summary and cleanup\n" | tee -a "$log"

	# identify how many wells pass both quality filters
	printf "well\tRNAseq_status\tWGBS_status\toverall_status\n" > "$outdir/summaries/all_statuses.txt" | tee -a "$log"	
	join -j 1 -o 1.1,2.1,1.2,2.2 -t $'\t' <( cut -f1,6 "$outdir/summaries/RNAseq_counts_summary.txt" | sort -k1,1 ) <( cut -f1,8 "$outdir/summaries/WGBS_methylation_summary.txt" | sort -k1,1 ) | awk -F$'\t' '{OFS=FS} {if ($3=="pass" && $4=="pass") {print $0,"pass"} else {print $0,"fail"}}' | cut -f2- >> "$outdir/summaries/all_statuses.txt"
	
	passQC=$( cut -f4 "$outdir/summaries/all_statuses.txt" | awk 'BEGIN{tt=0} $1=="pass" {tt+=1} END {print tt}' )
	echo " - A total of $passQC wells out of "$( echo $welllist | wc -w )" passed QC metrics for both RNA-seq and WGBS." | tee -a "$log"	
	
	echo " - Subsetting count matrices and methylation tiling matrices to keep only these wells..." | tee -a "$log"	
	mycols=$( awk -F$'\t' '{OFS=FS} $4 == "pass" {print $1}' "$outdir/summaries/all_statuses.txt" )
	cut -f$( echo "1-3,"$( grep -Fxn "$mycols" <( head -1 "$outdir/RNAseq_count_matrices/RNAseq_counts_QCpass.txt" | tr "\t" "\n" ) | cut -f1 -d":" | paste -sd, ) ) "$outdir/RNAseq_count_matrices/RNAseq_counts_QCpass.txt" > "$outdir/RNAseq_count_matrices/RNAseq_counts_QCpass_both.txt"
	cut -f$( echo "1-3,"$( grep -Fxn "$mycols" <( head -1 "$outdir/WGBS_tiling_matrices/WGBS_CpGme_over_${binsize}bp_bins_QCpass.txt" | tr "\t" "\n" ) | cut -f1 -d":" | paste -sd, ) ) "$outdir/WGBS_tiling_matrices/WGBS_CpGme_over_${binsize}bp_bins_QCpass.txt" > "$outdir/WGBS_tiling_matrices/WGBS_CpGme_over_${binsize}bp_bins_QCpass_both.txt"
	cut -f$( echo "1-3,"$( grep -Fxn "$mycols" <( head -1 "$outdir/WGBS_tiling_matrices/WGBS_CHGme_over_${binsize}bp_bins_QCpass.txt" | tr "\t" "\n" ) | cut -f1 -d":" | paste -sd, ) ) "$outdir/WGBS_tiling_matrices/WGBS_CHGme_over_${binsize}bp_bins_QCpass.txt" > "$outdir/WGBS_tiling_matrices/WGBS_CHGme_over_${binsize}bp_bins_QCpass_both.txt"
	cut -f$( echo "1-3,"$( grep -Fxn "$mycols" <( head -1 "$outdir/WGBS_tiling_matrices/WGBS_CHHme_over_${binsize}bp_bins_QCpass.txt" | tr "\t" "\n" ) | cut -f1 -d":" | paste -sd, ) ) "$outdir/WGBS_tiling_matrices/WGBS_CHHme_over_${binsize}bp_bins_QCpass.txt" > "$outdir/WGBS_tiling_matrices/WGBS_CHHme_over_${binsize}bp_bins_QCpass_both.txt"
	
	echo " - Cleaning up and deleting temp files..." | tee -a "$log"	
	# move all dedup'd BAM files to $outdir
	mkdir -p "$outdir/BAMfiles/RNAseq"; mkdir -p "$outdir/BAMfiles/WGBS"
	mv "$path_to_scratch/STAR/"*"_dedup"* "$outdir/BAMfiles/RNAseq"

	# compress all the bismark SAM files before moving
	numwells=$( ls -1 "$path_to_scratch/bismark/"*"_WGBS_dedup.sam" | wc -l )
	wellsperthread=$(( $(( numwells / threads )) + 1 ))
	ls -1 "$path_to_scratch/bismark/"*"_WGBS_dedup.sam" | split -a 2 -d -l $wellsperthread - "$path_to_scratch/bismarklist"
	
	pids=(); cmds=()
	for ff in "$path_to_scratch/bismarklist"*; do
		filelist=$( cat "$ff" | tr '\n' ',' )
		cmd="compress_sam_multi $filelist true" 
		eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
	done
	
	# wait for all processes to finish
	for ((i=0;i<${#pids[@]};++i)); do
		wait "${pids[i]}" || err_msg "SAM file sorting failed, command that failed: ${cmds[i]}" "$log"
	done
	
	# clean up temp files
	rm "$path_to_scratch/bismarklist"*
	mv "$path_to_scratch/bismark/"*"_dedup"* "$outdir/BAMfiles/WGBS"
	
	# also move over the per-position methylation BED files
	for context in CpG CHG CHH; do
		mkdir -p "$outdir/WGBS_per_position/$context"
		mv "$path_to_scratch/methylation_extractor/"*/*_"${context}.bed" "$outdir/WGBS_per_position/$context"
	done
	
	# delete all other temp files...
	[ "$keeptemp" = "false" ] && rm -rf "$path_to_scratch"

else
	printf "\nSkipping step 8. Wait, there are only 8 steps. You skipped the whole thing!" | tee -a "$log"
fi


te=$(date +%s); echo "Done at $(date). Time elapsed: $( displaytime $(($te - $ts)) )" | tee -a "$log"	

time_end=$(date)	# time run was started
time_es=$(date +%s)	# time run was started
echo "" | tee -a "$log"
echo "Run ended $time_end" | tee -a "$log"
echo "Total time elapsed: $( displaytime $(($time_es - $time_ss)) )" | tee -a "$log"


