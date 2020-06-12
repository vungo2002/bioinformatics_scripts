# this script will add a number of duplicates into a bam file by copying the reads and marking them as duplicates

from sys import argv
import sys 
import os
from collections import defaultdict
import pysam
import numpy as np

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]


### MAIN ###
inbam = argv[1]
duprate = float(argv[2])
outbam = argv[3]

if duprate > 0.9: 
    print("ERROR: duplicate rate is too high")
    sys.exit()
if duprate<0.0:
    print("ERROR: duplicate rate is negative")
    sys.exit()

mean_rate = duprate/(1-duprate)
r = 1 # first failure
dup_prob = mean_rate/(mean_rate + r)


infile = "../data/August_Illumina_WGS/1X_BAM_wo_dups_MPnormed/15-P4_S1.nodup.normed.1X.bam"

samfile = pysam.AlignmentFile(inbam, "rb")
pairedreads = pysam.AlignmentFile(outbam, "wb", template=samfile)

for read1, read2 in read_pair_generator(samfile):
    pairedreads.write(read1)
    pairedreads.write(read2)
    num_dup = np.random.negative_binomial(n=r, p=1-dup_prob) # determine the number of dups for this pair
    
    for i in range(num_dup):
        dupread1 = read1
        dupread2 = read2
        if dupread1.flag< 1024:
            dupread1.flag += 1024 # mark dup
        if dupread2.flag< 1024:
            dupread2.flag += 1024
        dupread1.query_name = read1.query_name + "_d"+str(i+1)
        dupread2.query_name = read2.query_name + "_d"+str(i+1)

        pairedreads.write(dupread1)
        pairedreads.write(dupread2)

pairedreads.close()
samfile.close()