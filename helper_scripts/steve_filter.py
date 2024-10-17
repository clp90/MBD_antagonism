#!/usr/bin/env python3

''' 
-------------------------
Version history:
	v.1.0: initial build	- 08/19/2020 by Colette L Picard
	v.2.0: 09/20/2020 by CLP
		- fixed Steve filter to only remove reads with 3 consecutive
		mCHHs -and- all positions between them also methylated, regardless of
		context (this is in line with Shawn Cokus' initial filter in his 2008 paper)
		- script now outputs reads that failed in a separate file
	v.3.0: 03/09/2021 by CLP
		- script can now run substantially faster if it's given a name-sorted input
		file (runs through file once instead of twice, which is default); must provide option
		--namesort to get runtime speed up (only applies to PE reads)
		- also added option to specify number of consecutive unmethylated CHHs to use for censoring
		(default is still 3, as before)
-------------------------
'''
 
import sys, os, argparse, re

if len(sys.argv) == 1:
	print("-------------------------")
	print("steve_filter v3.0		by Colette L. Picard, 03/09/2021")
	print("-------------------------")
	print("""Usage: steve_filter.py infile.sam outprefix
	
A very simple implementation of the Steve filter, which filters out bisulfite-sequencing
reads with consecutive unconverted Cs, with at least 3 of these in CHH context (since CHH 
methylation is rare, this is unlikely, and probably means the DNA wasn't properly bisulfite converted).

This script can handle single-end or paired-end data. Input file must be SAM, but doesn't
need to be sorted; however, script will run much faster with an input file sorted by read names
(set the --namesort option). Order of reads in input file will be preserved in output file. For PE data, 
if either mate of a pair is censored, its mate will also be removed. However, both
mates are evaluated separately (e.g. only censored if 3 consecutive mCHHs within that single read,
not across both mates).

Important note: this script recognizes paired-end mates by assuming they have identical read
names (first field in SAM file). If this is not the case, script will run, but paired-end reads 
will not be treated properly as paired!

Note for consecutive CHHs: the three CHHs must be methylated and consecutive but can have methylated
CGs or CHGs (but not unmethylated ones) between them too. For example: HHXH will be censored, but HHxH will not.

""")
	print("-------------------------")
	sys.exit(0)
	

# read in arguments
parser = argparse.ArgumentParser()
parser.add_argument('infile', help = 'SAM file of alignments')
parser.add_argument('outprefix', help = 'Name for output BAM file')
parser.add_argument('--namesort', default = False, action = "store_true",
	help = 'Input file has entries sorted by read name (speeds up script for PE data)')
parser.add_argument('--numCHHs', default = 3, type = int,
	help = 'Censor reads with this many consecutive unmethylated CHHs (default 3)')

args = parser.parse_args()

# process the numCHHs string
if args.numCHHs < 1: 
	print('Error, --numCHHs must be at least 1')
	sys.exit(1)	
if args.numCHHs > 10:
	print('Error, --numCHHs should be less than 10 (with such large values, no reads will get censored anyway)')
	sys.exit(1)	

regexstr = "H[XZ]*?" * (args.numCHHs - 1)+"H"

# -----------FUNCTIONS-------------
# function that takes a read and returns whether or not the read should be censored
def censorme(llstrip, methylidx):

	# find the 'methylation string' for this alignment
	if llstrip[methylidx].startswith('XM:Z:'):
		mestr = llstrip[methylidx].replace('.','')[5:]
	else:
		idxs = [i for i in llstrip if i.startswith('XM:Z:')] 
		if len(idxs) > 1:
			print('Strange error: more than one string starting with "XM:Z:" discovered in this alignment:')
			print(line)
			sys.exit(1)
		if len(idxs) == 0:
			print('Error: alignment',totreads,'did not contain the "XM" field, make sure your reads were aligned with bismark.')
			sys.exit(1)
		
		methylidx = llstrip.index(idxs[0])
		mestr = llstrip[methylidx].replace('.','')[5:]
				
	# should this read be censored according to the Steve filter?
	if bool(re.search(regexstr,mestr)):
		return True,methylidx
	else:
		return False,methylidx
		


# -----------MAIN-------------

# store number of reads that have been censored, and total reads processed
# (for reads mapped as pairs, only count first read of pair)

# Open SAM file and output file	
try:
	ff = open(args.infile, 'r')
except IOError as e:
	print('Could not open SAM file',args.infile)
	sys.exit(2)

try:
	oo = open(args.outprefix+'_pass.sam', 'w')
	of = open(args.outprefix+'_fail.sam', 'w')
except IOError as e:
	print('Could not open output file',args.outprefix,'for writing')
	sys.exit(2)

# skip header lines + add them to output file
line = ff.readline()
while line.startswith("@"):
	oo.write(line)
	of.write(line)
	line = ff.readline()

stevefilt = 0; totreads = 0; usecache = False; censorlist = {}; linenum = 0; methylstridx = 0

while line:
	linenum += 1
	if linenum % 100000 == 1:
		print('Reading line',linenum,'of SAM file...')
	ll = line.strip().split('\t')

	# check the SAM flag - is this alignment mapped? forward or reverse?
	flag = bin(int(ll[1]))[2:]; adjflag = '0'*(max(9-len(flag),0))+flag

	if adjflag[-3] == '1':							# read is mapped
		print('Error: SAM file contains unmapped entries (SAM flag is',ll[1],')')
		sys.exit(1)
		
	# is this read part of a pair?
	if adjflag[-1] == '1':
		# read is paired; if --namesort then mate should be next line, else cache
		if args.namesort == True:
			mate = ff.readline()
			matell = mate.strip().split('\t')
			
			if ll[0] != matell[0]:
				print('Error: --namesort option provided but mate of read:')
				print(ll[0])
				print('not found on next line. Next line has read:')
				print(matell[0])
				sys.exit(1)
						
			totreads += 1	
			res1,methylstridx = censorme(ll, methylstridx)
			res2,methylstridx = censorme(matell, methylstridx)
			
			if res1 == True or res2 == True:
				stevefilt += 1
				of.write(line)
				of.write(mate)
			else:
				oo.write(line)
				oo.write(mate)
		else:
			# reads are paired-end but not name-sorted, use cache and 2x pass approach instead
			usecache = True
			res1,methylstridx = censorme(ll, methylstridx)
			if res1 == True:
				censorlist[ll[0]] = 1
			
			# count only if this read is the first in the pair
			if adjflag[-7] == '1':
				totreads += 1
	else:
		# read is single-end, no need to cache or find mate so just output if ok
		totreads += 1
		res1,methylstridx = censorme(ll, methylstridx)
		if res1 == True:
			stevefilt += 1
			of.write(line)
		else:
			oo.write(line)

	line = ff.readline()


# If needed, do second pass, drop all lines in the censorlist
if usecache == True:
	ff.seek(0)
	print('Outputting all alignments that passed Steve filter (sort input file by name and use --namesort to speed this up)...')
	line = ff.readline()
	while line.startswith("@"):
		line = ff.readline()

	stevefilt = 0
	while line:
		ll = line.strip().split('\t')	
		if not ll[0] in censorlist:
			oo.write(line)
		else:
			of.write(line)
			stevefilt += 1
		line = ff.readline()
	
	stevefilt = stevefilt / 2

ff.close()
oo.close()
of.close()

# done, print summary
print("Done. A total of",totreads,"reads were processed (paired-end reads only count as one read)")
print(" - Of these,",stevefilt,"were censored by the --steve filter, and omitted from output file")
print(" -",totreads - stevefilt,"alignments remain after steve filter")



