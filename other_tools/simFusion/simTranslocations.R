#!/usr/bin/env Rscript
## ---------------------------
##
## Script name: simTranslocations.R
##
## Purpose of script: Use RSVSim to insert simulated fusions into wildtype hg19 reference genome 
##
## Author: Ankit Jambusaria
##
## Date Created: 10-7-21
##
## Email: ajambusaria@ambrygen.com
##
## ---------------------------
##
## Inputs: 
## 1. geneA: chr, start, end
## 2. geneB: chr, start, end 
##   
## Outputs:
## 1. genome_rearraged.fasta
## 2. translocations.csv
##
## ---------------------------

## Load Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RSVSim")
BiocManager::install("optparse")

library(optparse)
library(RSVSim)

## Read in command line args
option_list = list(
    ## geneA
    make_option(c("-a", "--geneAchr"), type="character", default=NULL, 
              help="Chromosome of geneA", metavar="character"),
    make_option(c("-b", "--geneAstart"), type="numeric", default=NULL, 
              help="gene A hg19 start coordinate", metavar="numeric"),
    make_option(c("-c", "--geneAend"), type="numeric", default=NULL, 
              help="gene A hg19 end coordinate", metavar="numeric"),
    ## geneB
    make_option(c("-d", "--geneBchr"), type="character", default=NULL, 
              help="Chromosome of geneB", metavar="character"),
    make_option(c("-e", "--geneBstart"), type="numeric", default=NULL, 
              help="gene B hg19 start coordinate]", metavar="numeric"),
    make_option(c("-f", "--geneBend"), type="numeric", default=NULL, 
              help="gene B hg19 end coordinate", metavar="numeric"),
    make_option(c("-n", "--fusionName"), type="character", default=NULL, 
              help="Fusion Name (e.g EML4-ALK)", metavar="character"),
    make_option(c("-o", "--outputDir"), type="character", default=NULL,
                help="Output Directory", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$geneAchr) | 
    is.null(opt$geneAstart) | 
    is.null(opt$geneAend) |
    is.null(opt$geneBchr) | 
    is.null(opt$geneBstart) | 
    is.null(opt$geneBend) |
    is.null(opt$fusionName)
   ){
    print_help(opt_parser)
    stop("All arguments must be supplied", call.=FALSE)
}

##
## ---------------------------
##
## Run RSVSim to generate fusion positive ref fasta
trans = GRanges(IRanges(opt$geneAstart,opt$geneAend), seqnames=opt$geneAchr,
chrB=opt$geneBchr, startB=opt$geneBstart, endB=opt$geneBend, balanced=TRUE)
names(trans) = opt$fusionName

#generate fusion containing ref fasta
print("Generating fusion containing reference...")
sim = simulateSV(output=opt$outputDir, chrs = c(opt$geneAchr, opt$geneBchr), regionsTrans=trans, bpSeqSize=50, random=FALSE)
print("Simulation Complete!")


