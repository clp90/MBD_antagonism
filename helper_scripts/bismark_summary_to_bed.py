#!/usr/bin/env python3

''' 
-------------------------
This script converts the output from bismark's methylation extractor function
into a more useful .bed file format. Specifically, the output contains a .bed
file with information about each cytosine in a particular context - # methylated
reads, # unmethylated and % methylated. See example below.

Can input as many input files as you want, so long as they are all for the 
same context (e.g. all have either X, Z or H). The output will contain the pooled
results for all of the provided input files. The last file in the input string is
the output file.

Usage:
python bismark_summary_to_bed.py infile1.txt [infile2.txt,...,infileN.txt] outfile

Bismark methylation summary output format example:
read	strand	chr	pos	me_call
WIGTC-HISEQ2:2:2303:18131:2366#GCCAAT/1;1	+	Chr4	4079224	Z
WIGTC-HISEQ2:2:2303:18131:2366#GCCAAT/1;1	+	Chr4	4079339	Z
WIGTC-HISEQ2:2:2303:3106:2367#GCCAAT/1;1	-	Chr5	18406211	z
WIGTC-HISEQ2:2:2303:3106:2367#GCCAAT/1;1	-	Chr5	18406219	z
WIGTC-HISEQ2:2:2303:3106:2367#GCCAAT/1;1	-	Chr5	18406226	z

Converted .bed file format (CpG file - have seperate file for each context):
Note .bed files use 0-based indexing while bismark output uses 1-based.
chr	start	end	#unme	#me	%me
chr4	4079223	4079224	0	1	100
chr4	4079338	4079339	0	1	100
chr5	18406210	18406211	1	0	0
chr5	18406218	18406219	1	0	0
chr5	18406225	18406225	1	0	0

Methylation call:
z	unmethylated CpG
Z	methylated CpG
x	unmethylated CHG
X	methylated CHG
h	unmethylated CHH
H	methylated CHH

v1.1	12/01/2014
by Colette Picard

Version history
v1.0 
- initial version (10/26/2014)
v1.1
- now works for all organisms instead of having the A. thaliana chromosomes hardcoded
- no longer outputs seperate files for ChrC, ChrM (can extract these manually w/ awk in later analyses)
- summaries of methylation levels on chrs are now unweighted averages over all sites w/ N > 5 reads
overlapping (previously used weighted mean over all sites on chr, leads to some sites disproportionately
affecting outcome)
- now outputs summary statistics (mean, sd) of read depth on each chromosome as well (in log file)
(uses Numpy v.1.9.0)
- reformatted log file into table form with the following fields (one record per chr):
	- chr = chromosome name
	- N_sites = total # of cytosines with at least 1 overlapping read
	- N_min5 = # of cytosines w/ at least 5 overlapping reads
	- mean_me_min5 = average methylation level across all sites with at least 5 reads (unweighted mean)
	- mean_depth = average depth (= # of reads overlapping a site)
	- min_depth = minimum depth (probably always 1)
	- max_depth = maximum depth observed
	- 25p_depth = depth at 25th percentile
	- median = depth at 50th percentile (median)
	- 75p_depth = depth at 75th percentile
v1.2
- checks that the methylation char is either z, x or h at the beginning instead of the end
- if none of the input files can be read, now says so instead of saying that " " is not z,x, or h

-------------------------
'''
 
import sys, numpy as np

if len(sys.argv)<=2: 
	if sys.argv[1] == "-h":
		print("Helpful information goes here.")
		sys.exit(1)
				
if len(sys.argv)<3: 
	print("Error: must specify at least one input file and output file.")
	sys.exit(1) 
	
outfile = sys.argv[len(sys.argv)-1]
infile_list = sys.argv[1:len(sys.argv)-1]

print("Running bismark_summary_to_bed v1.2 - by Colette L. Picard	06/16/2015")
print("Converting bismark output from the following files:")
for file in infile_list:
	print(file)
print("into a .bed file which will be saved as:")
print(outfile+".bed")

#-------------------------------------------------------------

char = ""	# store the methylation char for this set of input files
			# (that way, can check all inputs are, for example, for CHG context)

# For each input file, read the file & save the per chr, pos methylation info as a dict:
# me_scores[chr][pos] = [#me, #unme]
me_scores = {}
	
for infile in infile_list:
	print("\nProcessing input file",infile)
	# Read through the current input file
	try:
		f = open(infile, 'r') 
		line = f.readline()
		me_string = line.strip().split('\t')[4]		# get methylation string for this input file
		if char == "":
			print('Methylation string is',me_string)
			char = me_string
		elif char.lower() != "h" and char.lower() != "z" and char.lower() != "x":
			print('Error: methylation character',char,'is not z,x, or h!')
			print('Line is:',line)
			sys.exit(1)
		else:
			print('Checking if file has methylation string',char,'...')
			try:
				assert char.lower() == me_string.lower()
			except AssertionError:
				print('Error: input files do not have same methylation string. First file has',char,'while current file (',infile,') has',me_string)
				sys.exit(2)
			print('Passed! Current file has methylation string',char)
		
		while line:
			ll = line.strip().split('\t')	# split line by tabs
			assert len(ll) == 5				# should have 5 fields (see desc. above of inputs)
			# ll[2] = chr, ll[3] = pos, ll[4] = methylation call
			pos = int(ll[3])
			if not (ll[2] in me_scores):	# this chromosome/scaffold not yet encountered
				me_scores[ll[2]] = {}
		
			# store methylation info from this line
			if not (pos in me_scores[ll[2]]):
			# this position hasn't yet been put in dict, add it
				if str.isupper(ll[4]):
					me_scores[ll[2]][pos] = [1,0]	# methylated
				else:
					me_scores[ll[2]][pos] = [0,1]	# unmethylated
			else:
			# data already exists in dict for this position, add current record to it
				if str.isupper(ll[4]):
					me_scores[ll[2]][pos][0] += 1	# methylated
				else:
					me_scores[ll[2]][pos][1] += 1	# unmethylated
			
			line = f.readline()
		f.close()
	except IOError as e:
		print('Warning: could not open input file',infile,', skipping this file.')

if not me_scores:
	print("Error: no input files could be read")
	sys.exit(1)
print("\nDone processing input files. Outputting results to .bed file format...")

# print out to .bed file (unsorted - use sortBed to sort later) - also print
# summary statistics to log file
try:
	f_out = open(outfile+".bed",'w')
	log = open(outfile+"_log.txt",'w')
except IOError as e:
	print(e)
	print('Could not create output file',outfile,'.bed or log file',outfile,'_log.txt')
	sys.exit(2)

# output header line for log file
log.write("chr\tN_sites\tN_min5\tavg_methylation_min5\tmean_depth\tmin_depth\tmax_depth\t25p_depth\tmedian\t75p_depth\n")

# output main dataset to `outfile'.bed, also output summary stats
N_all=0; depth_all = []; f_min5_all = 0.0; N_min5_all = 0
for chr in me_scores:
	depth = []; f_min5 = 0.0; N = 0; N_min5 = 0
	for pos in me_scores[chr]:
		N+=1; depth.append(me_scores[chr][pos][0] + me_scores[chr][pos][1])		# add 1 to N sites, add read depth at that site to depth
		N_all+=1; depth_all.append(me_scores[chr][pos][0] + me_scores[chr][pos][1])
		if me_scores[chr][pos][0]+me_scores[chr][pos][1] >= 5:
			N_min5+=1; N_min5_all+=1
			f_min5 += float(me_scores[chr][pos][0])/(me_scores[chr][pos][0] + me_scores[chr][pos][1])		# if depth > 5, add frac methylated to running sum
			f_min5_all += float(me_scores[chr][pos][0])/(me_scores[chr][pos][0] + me_scores[chr][pos][1])	# if depth > 5, add frac methylated to running sum
		# write data to file:
		# chr	start	end	num_unmethylated	num_methylated	percent_methylated
		f_out.write(chr+'\t'+str(pos-1)+'\t'+str(pos)+'\t'+str(me_scores[chr][pos][1])+'\t'+str(me_scores[chr][pos][0])+'\t'+str(float(me_scores[chr][pos][0])/(me_scores[chr][pos][0]+me_scores[chr][pos][1])*100.0)+'\n')
	# output summary of methylation in current chr
	log.write(chr+'\t'+str(N)+'\t'+str(N_min5)+'\t')
	if N_min5 > 0:
		log.write(str(f_min5 / N_min5)+'\t')
	else:
		log.write('N/A\t')
	# also write summary stats for depth at each site
	d = np.array(depth)
	log.write(str(np.mean(d))+'\t'+str(np.min(d))+'\t'+str(np.max(d))+'\t')
	log.write(str(np.percentile(d,25))+'\t'+str(np.percentile(d,50))+'\t'+str(np.percentile(d,75))+'\n')

d = np.array(depth_all)		
if char.lower() == "z":
	log.write("overall_CpG\t")
elif char.lower() == "x":
	log.write("overall_CHG\t")
elif char.lower() == "h":
	log.write("overall_CHH\t")
else:
	print("Methylation char",char,"is not z,x, or h!")
	sys.exit(2)
	
if N_min5_all > 0:
	log.write(str(N_all)+'\t'+str(N_min5_all)+'\t'+str(f_min5_all / N_min5_all)+'\t'+str(np.mean(d))+'\t'+str(np.min(d))+'\t'+str(np.max(d))+'\t'+str(np.percentile(d,25))+'\t'+str(np.percentile(d,50))+'\t'+str(np.percentile(d,75))+'\n')
else:
	log.write(str(N_all)+'\t'+str(N_min5_all)+'\tN/A\t'+str(np.mean(d))+'\t'+str(np.min(d))+'\t'+str(np.max(d))+'\t'+str(np.percentile(d,25))+'\t'+str(np.percentile(d,50))+'\t'+str(np.percentile(d,75))+'\n')
	

f_out.close()
log.close()







