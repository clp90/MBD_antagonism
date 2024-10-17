#!/bin/bash

# Small function to pool N methylation datasets in BED-like format that were output from pseudobulk_WGBS.sh
# (which pseudobulks single cell DNA methylation data); output files from pseudobulk_WGBS.sh
# have files with seven fields: chr, start, end, unme, me, %me, # cells
# Adds an eighth column indicating number of plates (aka samples) that had reads overlapping that position
# Usage: pool_WGBS_bed.sh outprefix infiles
# Outputs a file called ${outprefix}_pool.bed

# NOTE: input files MUST BE SORTED!

pool_WGBS_bed ()
{
for ffid in `seq 2 $#`; do
echo "Processing ${!ffid}..."
if [ ! -f "${1}_pool.bed" ]; then
awk -F$'\t' '{OFS=FS} {print $0,1}' "${!ffid}" > "${1}_pool.bed"
else	
cat <( bedtools intersect -a "${1}_pool.bed" -b "${!ffid}" -wa -wb -sorted -f 1 -r | awk -F$'\t' '{OFS=FS} {print $1,$2,$3,$4+$12,$5+$13,($5+$13)/($4+$12+$5+$13)*100,$7+$15,$8+1}' ) <( bedtools intersect -a "${1}_pool.bed" -b "${!ffid}" -sorted -v ) <( bedtools intersect -b "${1}_pool.bed" -a "${!ffid}" -sorted -v | awk -F$'\t' '{OFS=FS} {print $0,1}' ) | sort -k1,1 -k2n,2 > "${1}_tmp.bed"
mv "${1}_tmp.bed" "${1}_pool.bed"
fi
done
}

[ -f "${1}_pool.bed" ] && { echo "Error: outfile ${1}_pool.bed already exists!"; exit 1; }

pool_WGBS_bed "$@"
