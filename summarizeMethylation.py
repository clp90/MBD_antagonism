#!/usr/bin/env python

'''
This function is called by compByFeature.sh - look at that script's description for a better
summary of the overall purpose of this pipeline and the role of this script within it.

Takes BED file with intersection of summarized methylation data (see scorePerPos/* in
Summarize_by_window.sh) with list of genomic features of interest, and outputs summarized
methylation within each feature. Specify which summarization metric should be used (default 1): 
	(1) Weighted mean methylation level (mean of %me at each site, weighted by seq depth at each site) 
	(2) Fraction of sites with methylation above a specified cutoff (default 0.5)

arguments:
	BED : BED file output from intersectBED using summarized methylation data and features of interest
	outfile : name of output file (including extension)
	method : (optional, default weighted mean) method to use to summarize methylation across the 
		regions of interest
	cutoff : (optional, only used if method == 2) cutoff to use for %methylated reads in order for
		site to be considered methylated (e.g. if == 0.5, at least 50% of reads at a site must be
		methylated for site as a whole to be considered methylated)
	minReads : (optional, default 0)

outputs:
	BED file containing chr, start, end, feature name (optional), score. 
	Feature name will only appear if in the input file, otherwise will
	contain '.' place marker.

usage:
	python summarizeMethylation.py BED.bed outfile [-method] [-cutoff] [-minReads]

'''
import sys, argparse, re

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument('BED', help = 'BED file containing intersected methylation and feature data.')
parser.add_argument('out', help = 'Output file name (with extension)')
parser.add_argument('-method', default = 1, type = int,
	help = 'Provide number of summarization metric you want to use: 1 = weighted mean, 2 = cutoff')
parser.add_argument('-cutoff', default = 0.5, type = float,
	help = 'If using cutoff (method == 2), give fraction methylation cutoff here (default 0.5)')
parser.add_argument('-minReads', default = 0, type = int,
	help = 'Sites with fewer than minReads reads are censored (sites == minReads kept).')
parser.add_argument('-omit', default = False, action = "store_true",
	help = 'Omit methylation value(s) in indicated interval of field 9 in merged file (must be in form e.g. Chr1:789-791)')
parser.add_argument('-qui', default = True, type = bool,
	help = 'Will suppress printing commentary to stdout - default.')
args = parser.parse_args()

if not args.qui:
	if args.method == 1:
		print('Summarizing methylation using weighted mean.')
	elif args.method == 2:
		print('Summarizing methylation by the fraction of sites with fraction methylated reads above',args.cutoff)

def importData(fname):
	'''
	inputs:
		fname - BED file with intersection of summarized methylation data (see scorePerPos/* in
		Summarize_by_window.sh) with genomic regions of interest.
	returns:
		meData - A dictionary of {feature's (chr, start, stop): {'name':feature name, 'sites': {start: (# un-me, # me)}}
	'''	
	# initialize the dictionary
	meData = {}
	# open and read file
	try:
		f = open(fname, 'r') 
	except IOError as e:
		print(e)
		print('Could not read input file',fname,'please check filename.')
		sys.exit(2)

	line = f.readline()
	# check if line is header, if so, skip
	if not 'Chr' and not 'chr' in line:
		line = f.readline()
	
	field9 = re.compile("^[A-Za-z]+[0-9]+:[0-9]+-[0-9]+$")
	
	while line:
		ll = line.strip().split('\t')	# split by tabs
		assert len(ll) >= 9				# must have at least 8 fields: site's chr, start, end, #me, #unme, frac
										# and feature's chr, start, end (optional name and strand)
		if int(ll[3])+int(ll[4]) >= args.minReads:
			# only keep sites with at least minReads reads total
			
			# if -omit switch on, test if field 9 in correct format
			proceed = True
			if args.omit:
				if field9.match(ll[9]):
					cchr = ll[9].split(':')[0]
					spos = ll[9].split(':')[1].split('-')[0]
					epos = ll[9].split(':')[1].split('-')[1]
			
					if ll[0] == cchr and ll[1] >= spos and ll[2] <= epos:		# this record is in the interval to skip
#						print "This record is in interval to skip; skipping:"
#						print ll
						proceed = False

			if proceed and (not (ll[0],int(ll[7]),int(ll[8])) in meData):
				# this genomic feature not yet in dictionary
				# check if name provided (== ll[9])
				if len(ll) >= 10:
					meData[(ll[0],int(ll[7]),int(ll[8]))] = {'name':ll[9],'sites':{int(ll[1]):(int(ll[3]),int(ll[4]))}}
				else:
					meData[(ll[0],int(ll[7]),int(ll[8]))] = {'name':'','sites':{int(ll[1]):(int(ll[3]),int(ll[4]))}}
			elif proceed:
				meData[(ll[0],int(ll[7]),int(ll[8]))]['sites'][int(ll[1])] = (int(ll[3]),int(ll[4]))
		
		line = f.readline()
			
	f.close()
	return meData
	
def getWeightedMean(meSites):
	'''
	Returns the weighted mean of the CpG, CHH or CHG sites in meSites.
	inputs:
		meSites - A dictionary of {start: (# un-me, # me)}
	returns:
		stat - single # summarizing the methylation status of this feature
	'''
	num_me = 0		# total number methylated reads
	num_unme = 0	# total number unmethylated reads
	for site in meSites:
		num_unme += meSites[site][0]
		num_me += meSites[site][1]
	
	return float(num_me)/float(num_me + num_unme)
	
def getFractionMethyl(meSites, cutoff):
	'''
	Returns the fraction of CpG, CHH or CHG sites in meSites that have > cutoff fraction
	of reads that are methylated (e.g. if cutoff = 0.5, returns the number of sites with at least
	half the reads showing methylation at that site).
	inputs:
		meSites - A dictionary of {start: (# un-me, # me)}
		cutoff - fraction methylation cutoff for site to be considered methylated (0 <= cutoff <= 1)
	returns:
		stat - fraction of sites (for which we have data) with methylation levels above the cutoff
	'''
	numMeSites = 0
	numUnmeSites = 0
	assert cutoff <= 1.0 and cutoff >= 0.0
	
	for site in meSites:
		if float(meSites[site][1])/float(meSites[site][0]+meSites[site][1]) >= cutoff:
			numMeSites+=1
		else:
			numUnmeSites+=1
			
	return float(numMeSites)/float(numMeSites+numUnmeSites)
	
	
def main():
	# import the data
	if not args.qui:
		print('Importing data from',args.BED)
	meData = importData(args.BED)
	
	# for each feature, calculate the desired methylation summary statistic and output to file
	if not args.qui:
		print('Calculating summarized methylation levels...')
	out = open(args.out,'w')
	for feature in meData:
		if args.method == 1:
			# weighted mean
			stat = getWeightedMean(meData[feature]['sites'])
		elif args.method == 2:
			# fraction above cutoff
			stat = getFractionMethyl(meData[feature]['sites'],args.cutoff)
		else:
			print('This summarizing method does not exist or is not implemented. Check available methods in file description.')
			sys.exit(2)
		
		# write output to file (chr, feature start, feature end, feature name, summary stat, N)
		# note: N = number of sites that can be methylated in this feature
		out.write(feature[0]+'\t'+str(feature[1])+'\t'+str(feature[2])+'\t'+meData[feature]['name']+
			'\t'+str(stat)+'\t'+str(len(meData[feature]['sites']))+'\n')
		
	out.close()
	if  not args.qui:	
		print('Done.')
	
main()	
	
	
	
	
