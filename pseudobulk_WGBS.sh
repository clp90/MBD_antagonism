#!/usr/bin/env bash
trap 'pkill -P $$; exit' SIGINT SIGTERM					# kills child processes on exit

# ------------------------------------------------------------------------------------
# v1.0 by Colette Picard
# 04/20/2022
# ------------------------------------------------------------------------------------
	
# Usage:
# pseudobulk_WGBS.sh [options] -d BEDdir -c clusters.txt -o outprefix

# -------------------------
# Version history:
# v.1.0: initial build - 04/20/2022 by CLP
# -------------------------

# To be added:
# - TBD

# Description printed when "help" option specified:
read -d '' usage <<"EOF"
v1.0 by Colette Picard, 04/20/2022
This is a companion script to the mapping script for snmCAT-seq (snmCT_seq_map.sh)
which can take the WGBS output of that script and, in combination with a set of provided
nuclei clusters, 'pseudobulk' the WGBS data by pooling all nuclei in each cluster.

Takes advantage of threads if available, should almost linearly speed up analysis.
User does not need to specify threads available. One thread is used per cluster/context.

Inputs:
- BEDdir should point to the WGBS_per_position/ folder in the snmCAT_seq_map.sh output dir, or more
generally to a folder containing methylation data in BED-like format (chr, start, end, unme, me, %me)
in subdirectories by context, with filenames *_{well name}_{context}.bed
- clusters should be a tab-separated file w/ 2 columns corresponding to well name and cluster ID, 
with no header. For example:
A1	SN
A3	unk_1
A7	veg
A11	mid_VN
A13	veg
- a chrom.sizes file of genome used for alignment, for example:
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
ChrC	154478
ChrM	366924
- outprefix is prefix for all output files and temp files

NOTE: assumes that 6th field of input BED files have values in [0,100] (e.g. percents, not fractions)
This is always true if you're using snmCAT_seq_map.sh output, but keep in mind. If this isn't true for input
files, you'll get some really odd results.
- Also, assumes all input files are sorted by position (also always true if you're using snmCAT_seq_map.sh output)


EXTRA STUFF:
-------------------------------
TBD
-------------------------------

Usage:
pseudobulk_WGBS.sh [options] -d BEDdir -c clusters.txt -o outprefix

User-specified options (defaults indicated in brackets):
Required arguments:
	-d BEDdir : location of the WGBS_per_position/ folder in the snmCAT_seq_map.sh output dir
	-c clusters : tab-separated file w/ 2 columns corresponding to well name and cluster ID
	-f chromsizes : chrom.sizes file of genome used to align reads
	-o outprefix : prefix for all output files and temp files
Additional options:
	-m minsites : minimum number of sites that must be covered for well to be used (across all contexts) [500]
	-x minreads : minimum total number of WGBS reads in library for well to be used (across all contexts) [5000]
	-t contexttouse : only analyze indicated context (default analyzes all contexts in file) []
Flag options:
	-r : allow overwrite of existing files (WARNING: if you have other files starting w/ $outprefix they may be deleted!) [overwrite=false]
	-0 : checks that all required programs installed on PATH and all required helper scripts can be located, then exits without running
	-h : prints this version and usage information

Must be installed on your PATH:
	- bedtools (v.2.30.0 by quinlanlab.org and others)
	- bedGraphToBigWig (part of UCSC tools, tested on v4)
	
------------------------------------------------------------------------------------
EOF

[[ $# -eq 0 ]] && { printf "%s\n" "$usage"; exit 0; } 		# if no user-supplied arguments, print usage and exit

# ----------------------
# Get user-specified arguments
# ----------------------

# Initiate environment
scriptDir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )	# location of this script 
workdir=$( pwd )											# working directory

# Additional options:
# ----------------------
minsites=500								# minimum number of sites that must be covered for well to be used (across all contexts)
minreads=5000								# minimum total number of WGBS reads in library for well to be used (across all contexts)
contexttouse=""							# only analyze indicated context (default analyzes all contexts in file)

# Flag options:
# ----------------------
overwrite=false							# allow overwrite of existing files (WARNING: if you have other files starting w/ $outprefix they may be deleted!)

# Required arguments:
# ----------------------
BEDdir=""							# location of the WGBS_per_position/ folder in the snmCAT_seq_map.sh output dir
clusters=""							# tab-separated file w/ 2 columns corresponding to well name and cluster ID
chromsizes=""							# chrom.sizes file of genome used to align reads
outprefix=""							# prefix for all output files and temp files

checkdep=false

# ----------------------
while getopts "d:c:f:o:m:x:t:r0h" opt; do
	case $opt in
		d)	# location of the WGBS_per_position/ folder in the snmCAT_seq_map.sh output dir
			BEDdir="$OPTARG"
			;;
		c)	# tab-separated file w/ 2 columns corresponding to well name and cluster ID
			clusters="$OPTARG"
			;;
		f)	# chrom.sizes file of genome used to align reads
			chromsizes="$OPTARG"
			;;
		o)	# prefix for all output files and temp files
			outprefix="$OPTARG"
			;;
		m)	# minimum number of sites that must be covered for well to be used (across all contexts)
			minsites="$OPTARG"
			;;
		x)	# minimum total number of WGBS reads in library for well to be used (across all contexts)
			minreads="$OPTARG"
			;;
		t)	# only analyze indicated context (default analyzes all contexts in file)
			contexttouse="$OPTARG"
			;;
		r)	# allow overwrite of existing files (WARNING: if you have other files starting w/ $outprefix they may be deleted!)
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
command -v bedtools >/dev/null 2>&1 || { echo "Error: bedtools is required on PATH but was not found"; exit 1; }
command -v bedGraphToBigWig >/dev/null 2>&1 || { echo "Error: bedGraphToBigWig is required on PATH but was not found"; exit 1; }

# Done checking all requirements. Stop here if -0 flagged.
# ----------------------
"$checkdep" && exit 0

# Check all required inputs are provided
# ----------------------
[ -z "$BEDdir" ] && { echo "Error: -d BEDdir is a required argument (location of the WGBS_per_position/ folder in the snmCAT_seq_map.sh output dir)"; exit 1; }
[ -z "$clusters" ] && { echo "Error: -c clusters is a required argument (tab-separated file w/ 2 columns corresponding to well name and cluster ID)"; exit 1; }
[ -z "$chromsizes" ] && { echo "Error: -f chromsizes is a required argument (chrom.sizes file of genome used to align reads)"; exit 1; }
[ -z "$outprefix" ] && { echo "Error: -o outprefix is a required argument (prefix for all output files and temp files)"; exit 1; }

# Check that all inputs can be located
# ----------------------
[ -d "$BEDdir" ] || { echo "Error: directory $BEDdir couldn't be located or isn't a directory"; exit 1; }
contexts=$( ls "$BEDdir" | tr '\n' ' ' )
contextlist=( $contexts )
if [ ! -z "$contexttouse" ]; then
	[[ ! " ${contexts[*]} " =~ " ${contexttouse} " ]] && { echo "Error: context provided to -t option (${contexttouse}) not found in $BEDdir (which contains contexts ${contexts[*]})"; exit 1; }
	contexts="$contexttouse"
	contextlist=( $contexts )
fi
numcontexts="${#contextlist[@]}"
[ -f "$chromsizes" ] || { echo "Error: could not open file $chromsizes"; exit 1; }
[ -f "$clusters" ] || { echo "Error: could not open file $clusters"; exit 1; }


# ----------------------
# Helper functions for this script
# ----------------------
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

pool_WGBS_bed ()
# Small function to pool N methylation datasets in BED format
# Adds a seventh column indicating number of wells that had reads overlapping that position
# Usage: pool_WGBS_bed outprefix infiles
{
	[ -f "${1}_pool.bed" ] && { echo "Internal error in pool_WGBS_bed(); outfile ${1}_pool.bed already exists!"; exit 1; }
	for ffid in `seq 2 $#`; do
		if [ ! -f "${1}_pool.bed" ]; then
			awk -F$'\t' '{OFS=FS} {print $0,1}' "${!ffid}" > "${1}_pool.bed"
		else	
			cat <( bedtools intersect -a "${1}_pool.bed" -b "${!ffid}" -wa -wb -f 1 -sorted -r | awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$4+$11,$5+$12,($5+$12)/($4+$11+$5+$12)*100,$7+1}' ) <( bedtools intersect -a "${1}_pool.bed" -b "${!ffid}" -v -sorted ) <( bedtools intersect -b "${1}_pool.bed" -a "${!ffid}" -v -sorted | awk -F$'\t' '{OFS=FS} {print $0,1}' ) | sort -k1,1 -k2n,2 > "${1}_tmp.bed"
			mv "${1}_tmp.bed" "${1}_pool.bed"
		fi
	done
}


# ----------------------
# Main code
# ----------------------
log="${outprefix}_log.txt" 	# create log file
time_start=$(date)				# time run was started
time_ss=$(date +%s)				# time run was started (in seconds)

# Output user-derived options to stdout and to log file
# ----------------------
echo "" | tee -a "$log"
echo "Running pseudobulk_WGBS v1.0 by Colette Picard (4/20/2022):" | tee -a "$log"
echo "Run start on: $time_start" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Working directory: $( pwd )" | tee -a "$log"
echo "WBGS data folder: $BEDdir" | tee -a "$log"
echo " - Detected $numcontexts contexts: $contexts" | tee -a "$log"
echo "Clusters file: $clusters" | tee -a "$log"
echo "chrom.sizes file: $chromsizes" | tee -a "$log"
echo "Prefix for output files: $outdir" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Additional settings:" | tee -a "$log"
echo "Nuclei filtering parameters:" | tee -a "$log"
echo " - Min number of sites covered: $minsites" | tee -a "$log"
echo " - Min aligned reads: $minreads" | tee -a "$log"
echo "-------------------------" | tee -a "$log"
echo "Full command used:" | tee -a "$log"
echo "$0 $@" | tee -a "$log"
echo "-------------------------" | tee -a "$log"


# ----------------------
# Step 1: read in clusters file and check all files can be located in BEDdir
# ----------------------
echo " - Checking inputs and identifying wells that pass QC for each cluster"

[ -f "$clusters" ] || { echo "Error: could not open clusters file $clusters"; exit 1; }
cluslist=( $( cut -f2 $clusters | sort | uniq -c | sed 's/^ *//g' | cut -f2 -d' ' ) )
origcounts=( $( cut -f2 $clusters | sort | uniq -c | sed 's/^ *//g' | cut -f1 -d' ' ) )

for ((i=0;i<${#cluslist[@]};++i)); do
	for context in $contexts; do
		if "$overwrite"; then
			rm -f "${outprefix}_${cluslist[i]}_${context}_filelist.txt"
		else
			[ -f "${outprefix}_${cluslist[i]}_${context}_filelist.txt" ] && { echo "Error: temp file ${outprefix}_${cluslist[i]}_${context}_filelist.txt already exists; either delete or allow automatic overwrite with the -r option"; exit 1; }
		fi
	done
done

[ -f "${outprefix}_per_well_stats.txt" ] && [ "$overwrite" = "false" ] && { echo "Error: output file ${outprefix}_per_well_stats.txt already exists; either delete or allow automatic overwrite with the -r option"; exit 1; }
printf "well\tcluster\ttotreads\ttotsites\tpassQC\n" > "${outprefix}_per_well_stats.txt"

while IFS= read -r line; do
	well=$( echo "$line" | cut -f1 )
	clus=$( echo "$line" | cut -f2 )
	[ "$(( $i % 20 ))" -eq 0 ] && echo "  - Processing ${i}th well"
	
	totreads=0; totsites=0; passQC=false
	for context in $contexts; do
		ff=$( ls "$BEDdir/$context/"*"_${well}_${context}.bed" )
		[ $( echo "$ff" | wc -w ) -gt 1 ] && { echo "Error: multiple files for ${well} detected in ${context}, check your inputs"; exit 1; }
		[ $( echo "$ff" | wc -w ) -eq 0 ] && { echo "Error: could not find $context methylation BED file for ${well} in $BEDdir"; exit 1; }
		
		# get status to determine if passes QC
		res=$( awk -F$'\t' '{unme+=$4; me+=$5; sites+=1} END {print unme,me,sites}' "$ff" )
		unme=$( echo "$res" | cut -f1 -d' ' ); me=$( echo "$res" | cut -f2 -d' ' ); ss=$( echo "$res" | cut -f3 -d' ' )
		trs=$(( unme + me ))
		
		totreads=$(( totreads + trs ))
		totsites=$(( totsites + ss ))
	done
	
	# passes QC?
	[ "$totreads" -ge "$minreads" ] && [ "$totsites" -ge "$minsites" ] && passQC=true
	
	# if passes QC, add to filelists	
	if "$passQC"; then
		for context in $contexts; do	
			ff=$( ls "$BEDdir/$context/"*"_${well}_${context}.bed" )	
			echo "$ff" >> "${outprefix}_${clus}_${context}_filelist.txt"
		done
	else
		echo "Well $well failed QC"
	fi
	
	# add to summary file
	printf "$well\t$clus\t$totreads\t$totsites\t$passQC\n" >> "${outprefix}_per_well_stats.txt"
	
	i=$(( i + 1 ))	
done < "$clusters"

# summarize filtering results
context1=$( echo "$contexts" | cut -f1 -d' ' )
echo " - Done, summary:"
echo "    cluster	tot_wells	passQC_wells" 
for ((i=0;i<${#cluslist[@]};++i)); do
	if [ -f "${outprefix}_${cluslist[i]}_${context1}_filelist.txt" ]; then	
		passqcnum=$( cat "${outprefix}_${cluslist[i]}_${context1}_filelist.txt" | wc -l )
		echo "    ${cluslist[i]}	${origcounts[i]}	$passqcnum"
	else
		echo "    ${cluslist[i]}	${origcounts[i]}	0"
	fi
done


# ----------------------
# Step 2: combine methylation data for each cluster
# ----------------------
echo " - Combining methylation data for each cluster"

pids=(); cmds=()
for ((i=0;i<${#cluslist[@]};++i)); do
	for context in $contexts; do
		oop="${outprefix}_${cluslist[i]}_${context}"
		if [ -f "${outprefix}_${cluslist[i]}_${context}_filelist.txt" ]; then		
			if "$overwrite"; then
				rm -f "${oop}_pool.bed"
			else
				[ -f "${oop}_pool.bed" ] && { echo "Error: output file ${oop}_pool.bed already exists; either delete or allow automatic overwrite with the -r option"; exit 1; }
			fi
			cmd="pool_WGBS_bed $oop $( cat ${outprefix}_${cluslist[i]}_${context}_filelist.txt | tr '\n' ' ' )"
#			echo "Starting process $cmd"
			eval "$cmd" & pids+=( $! ); cmds+=( "$cmd" )
		else
			echo "Skipping ${cluslist[i]} ${context} which had no wells passing QC"
		fi
	done
done

# wait for all processes to finish
for ((j=0;j<${#pids[@]};++j)); do
#	echo "Waiting for process ${pids[j]}"
	wait "${pids[j]}" && echo "Process ${pids[j]} finished successfully" || { echo "pooling WGBS data failed, command that failed: ${cmds[j]}"; exit 1; }
done

# delete filelist temp files
for ((i=0;i<${#cluslist[@]};++i)); do
	for context in $contexts; do
		[ -f "${outprefix}_${cluslist[i]}_${context}_filelist.txt" ] && rm "${outprefix}_${cluslist[i]}_${context}_filelist.txt"
	done
done


# ----------------------
# Step 3: convert methylation files to .bw files and summarize
# ----------------------
echo " - Converting to bigwig"

pids=(); cmds=()
for context in $contexts; do
	for ((i=0;i<${#cluslist[@]};++i)); do
		oop="${outprefix}_${cluslist[i]}_${context}"
		if [ -f "${oop}_pool.bed" ]; then
			# for each, convert to bedGraph; make a version with zeros replaced by -10; then convert both to bigWig
			# make these as temp scripts then run each as subprocess
			echo '#!/usr/bin/env bash' > "${oop}_tmp.sh"
			chmod 755 "${oop}_tmp.sh"
		
			echo "cut -f1-3,6 \"${oop}_pool.bed\" > \"${oop}_pool.bedgraph\"" >> "${oop}_tmp.sh"
			echo "awk -F\$'\t' '{OFS=FS} \$4==0{\$4=-10}1' \"${oop}_pool.bedgraph\" > \"${oop}_pool_forIGV.bedgraph\"" >> "${oop}_tmp.sh"
			echo "bedGraphToBigWig \"${oop}_pool.bedgraph\" \"${chromsizes}\" \"${oop}_pool.bw\"" >> "${oop}_tmp.sh"
			echo "bedGraphToBigWig \"${oop}_pool_forIGV.bedgraph\" \"${chromsizes}\" \"${oop}_pool_forIGV.bw\"" >> "${oop}_tmp.sh"
			echo "" >> "${oop}_tmp.sh"
		
			"./${oop}_tmp.sh" & pids+=( $! ); cmds+=( "${oop}_tmp.sh" )	
		else
			echo "${oop}_pool.bed not found, skipping"
		fi
	done
done	
# wait for all processes to finish
for ((j=0;j<${#pids[@]};++j)); do
	wait "${pids[j]}" && echo "Process ${pids[j]} finished successfully" || { echo "bigwig conversion failed, command that failed: ${cmds[j]}"; exit 1; }
done

# delete the .bedgraph and .sh files
for ((i=0;i<${#cluslist[@]};++i)); do
	for context in $contexts; do
		oop="${outprefix}_${cluslist[i]}_${context}"
		if [ -f "${oop}_pool.bed" ]; then
			rm "${oop}_tmp.sh" "${oop}_pool.bedgraph" "${oop}_pool_forIGV.bedgraph"
		fi
	done
done

# summarize average methylation per chromosome in each cluster
echo " - Summarizing methylation levels per cluster"
[ -f "${outprefix}_per_cluster_stats.txt" ] && [ "$overwrite" = "false" ] && { echo "Error: output file ${outprefix}_per_cluster_stats.txt already exists; either delete or allow automatic overwrite with the -r option"; exit 1; }
printf "cluster\tcontext\tchr\ttot_unme\ttot_me\tfrac_me\tsites_cov\tavg_cells_per_site\n" > "${outprefix}_per_cluster_stats.txt"
chrlist=$( cut -f1 "$chromsizes" | tr '\n' ' ' )
for ((i=0;i<${#cluslist[@]};++i)); do
	echo "   - Summarizing ${cluslist[i]}..."
	for context in $contexts; do
		oop="${outprefix}_${cluslist[i]}_${context}"
		if [ -f "${oop}_pool.bed" ]; then
			for chr in $chrlist; do
				res=$( awk -F$'\t' -v cc="$chr" '{OFS=FS} $1==cc {unme+=$4; me+=$5; cells+=$7; sites+=1} END {print cc,unme,me,me/(unme+me),sites,cells/sites}' "${oop}_pool.bed" )
				printf "${cluslist[i]}\t${context}\t$res\n" >> "${outprefix}_per_cluster_stats.txt"
			done
		fi
	done
done

time_end=$(date)	# time run was started
time_es=$(date +%s)	# time run was started
echo "" | tee -a "$log"
echo "Run ended $time_end" | tee -a "$log"
echo "Total time elapsed: $( displaytime $(($time_es - $time_ss)) )" | tee -a "$log"













































