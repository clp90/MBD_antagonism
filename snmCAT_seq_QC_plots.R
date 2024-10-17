#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(reshape2))

# version 1.0 (06/21/2022)
# -------------------------
# Version history:
# v.1.0: initial build - 06/21/2022
# -------------------------

# Description: helper script for snmCAT_seq_map.sh that generates QC plots. Current supported modes (set by --mode):
# (1) plate quality plot - shows a representation of the 384 well plate, labelling wells that passed QC or failed QC
#		for either RNA-seq or WGBS



# Usage: boxplot.R [options] infile.txt outfile.pdf
parser = ArgumentParser()
parser$add_argument("infile", nargs=1, help = "tab-delimited input file, see description for details")
parser$add_argument("outfile", nargs=1, help = "name for output plot (in PDF format)")
parser$add_argument("--mode", type="integer", default=1, help = "mode to use/type of plot to make, see usage info (call without arguments), default 1 (plate quality plot)")

args <- commandArgs(trailingOnly = TRUE)	
if (length(args) <= 1) {
	cat("Usage: boxplot.R [options] infile.txt outfile.png\n")
	cat("----------------------\n")
	cat("infile.txt is the tab-delimited input file, see description below for details\n")
	cat("outfile.png is the name for output plot (in PDF format)\n")
	cat("----------------------\n")
	cat("Options:\n")
	cat("--mode <int> : mode to use/type of plot to make, see usage info (call without arguments), default 1 (plate quality plot)\n")
	cat("----------------------\n")
	cat("\n")
	cat("Description:\n")
	cat("Makes QC plots as a helper script for snmCAT_seq_map.sh. Currently supported plots (set by --mode) and example usage:\n")
	cat("\n")
	cat(" (1) Plate quality plot: shows which wells in the plate passed QC and which failed, useful for seeing if there\n")
	cat("     are any issues in library prep causing certain wells to fail consistently. Expected input is the summaries/all_statuses.txt\n")
	cat("     file output by snmCAT_seq_map.sh, with this format:\n")
	cat("well	RNAseq_status	WGBS_status	overall_status\n")
	cat("A1	pass	fail	fail\n")
	cat("A10	pass	pass	pass\n")
	cat("etc.\n")
	cat("\n")
	cat("----------------------\n")
}

opt <- parser$parse_args()

infile = opt$infile
outfile = opt$outfile

if (opt$mode == 1) {
	alldata = read.table(infile, header=TRUE, sep="\t", stringsAsFactors = FALSE)
	alldata$col = as.numeric(substr(alldata$well, 2, length(alldata$well)))
	alldata$row = substr(alldata$well, 1, 1)
	ll <- letters[1:26]
	alldata$rownum = 17 - as.numeric(match(tolower(alldata$row), ll))
	
	alldata$colortouse = ifelse(alldata$overall_status == "pass", "chartreuse", ifelse(alldata$RNAseq_status == "pass", "darkgoldenrod", ifelse(alldata$WGBS_status == "pass", "darkorange", "red")))
	
	pdf(outfile, width = 9, height = 5)
	p = ggplot(alldata, aes(x = col, y = rownum)) + geom_point(color = alldata$colortouse, size = 8) + theme_bw() + xlab("") + ylab("") + scale_x_continuous(labels = as.character(alldata$col), breaks = alldata$col) + scale_y_continuous(labels = as.character(alldata$row), breaks = alldata$rownum)


	print(p)
	graphics.off()


}

