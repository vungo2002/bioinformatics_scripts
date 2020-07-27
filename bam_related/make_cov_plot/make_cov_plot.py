#!/usr/bin/env python3

## TODO: check to see if bam file is indexed 

"""
Input: a bam file, option to ignore duplicates or not 
Output: a coverage plot of said bam file, over all chromosomal regions (hg38)
"""

import os
from sys import argv
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import numpy as np

bamfile = argv[1]
ignoredup = argv[2] # option to ignore duplicates or not, two values --nodup --withdup
outfile = argv[3] # a plot in pdf format 

if ignoredup != "--nodup" and ignoredup != "--withdup":
    print("Please choose --nodup or --withdup")
    sys.exit(-1)
if ignoredup == "--nodup":
    ignoredup = True
    print("Plotting Without Duplicates")
else:
    ignoredup = False
    print("Plotting WITH Duplicates")

# first, generate the coverage bedgraph from the bam file (using deeptools)
# if the bedgraph file exists, skip this step
tmp = bamfile.split("/")[:-1]
if len(tmp) < 1:
    wd = "./"
else:
    wd = '/'.join(tmp) + '/'
print(wd)
#sys.exit()


covfile = wd + bamfile.split("/")[-1].replace(".bam",".cov.bedgraph")
if os.path.exists(covfile):
    # already got the cov file, skip the deeptools part
    print("found already calculated coverages :)")
    pass
else:
    print("Calculating the coverages using deeptools bamCoverage")
    if ignoredup:
        cmd = "bamCoverage --bam %s -o %s --binSize 50000 --outFileFormat=bedgraph --numberOfProcessors=2 --ignoreDuplicates" %(bamfile, covfile)
    else:
        cmd = "bamCoverage --bam %s -o %s --binSize 50000 --outFileFormat=bedgraph --numberOfProcessors=2" %(bamfile, covfile)
    os.system(cmd)




# then read the bedgraph file and produce the plot 
infile = covfile
coverages = []
x_labels = []
for line in open(infile):
    tmp = line.strip().split()
    cov = float(tmp[-1])
    chromname = tmp[0]
    if "_" in chromname or "chr" not in chromname or "chrM" in chromname or "EBV" in chromname:
        continue
    #if tmp[0] != "chr1":
    #    continue
    coverages += [cov]
    x_labels += [tmp[0]]
len(coverages)
# norm it down to larger bins:
merge_factor = 5 # merge the bins into even larger bins
readsize = 150
covs = []
curtotal = 0
x_labels_new = []

for i in range(len(coverages)):
    curtotal += coverages[i]
    i += 1
    if i%merge_factor == 0: # average per 100 bins -> total length = 100*binsize = 100*5k = 500k
        covs += [curtotal/merge_factor/50000*readsize]
        curtotal = 0
        x_labels_new += [x_labels[i-1]]


# normalize the coverage by mean:
average_cov = np.mean(covs)
for i in range(len(covs)):
    covs[i] = covs[i]/average_cov*2 # normalize to 2n 


x_labels_locs = {}
x_labels_avg_covs = {}
chromorder = []
for i in range(len(x_labels_new)):
    chrom = x_labels_new[i]
    if "chr" not in chrom or "_" in chrom: # ignore the junks
        continue
    if chrom not in x_labels_locs:
        x_labels_locs[chrom] = 1
        chromorder += [chrom]
    else:
        x_labels_locs[chrom] += 1

# create locations for the ticks
total_counts = 0
for i in x_labels_locs:
    total_counts += x_labels_locs[i]
total_counts
newlocs = []
prev_loc = 0
boundary_locs = []
for chrom in chromorder:
    chromsize = x_labels_locs[chrom]
    newlocs += [prev_loc + chromsize/2]
    boundary_locs += [prev_loc + chromsize]
    prev_loc = boundary_locs[-1]
for i in range(len(newlocs)):
    print (newlocs[i], chromorder[i])




# ploting 
print("plotting", outfile,"...")
fig, ax = plt.subplots(nrows=1, ncols=1 )
ax.scatter(x=range(len(covs)), y = covs, alpha=0.5, color = "green", s = 0.5)
fig.set_figheight(10)
fig.set_figwidth(20)
ax.set_xticks(ticks=newlocs)
ax.set_xticklabels(labels=chromorder, rotation=90)
ax.set_title(bamfile.split('/')[-1])

plt.ylim(0,10)
for loc in boundary_locs:
    plt.axvline(x=loc, ms = 0.2, alpha=0.7, color='black')

plt.savefig(outfile)
#plt.show()