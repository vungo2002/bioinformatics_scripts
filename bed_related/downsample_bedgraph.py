#!/usr/bin/env python3 

from sys import argv
import sys
from scipy.stats import binom

"""
This script takes:
- An input bedgraph file with with 4th column being the number of reads in the region
- A desired total number of reads 
- Output bedgraph file
"""

in_bg = argv[1]
desired_total_num = int(argv[2])
out_bg = argv[3]


# calculate total number of reads in the input file
read_vector = []
for line in open(in_bg):
    line = line.strip().split('\t')
    read_vector += [int(line[3])]
in_total_reads = sum(read_vector)

if desired_total_num > in_total_reads:
    print("Error, total number of reads is {}, which is lower than the desired number {}".format(str(in_total_reads), 
        str(desired_total_num)))
    sys.exit()

ds_factor = desired_total_num/in_total_reads
print("Down-sampling factor is {}".format(str(ds_factor)))

total_out_read = 0

with open(out_bg, 'w') as f:
    for line in open(in_bg):
        line = line.strip().split("\t")
        cov = int(line[3])
        new_read_cov = binom.rvs(cov, ds_factor)
        total_out_read += new_read_cov
        line[3] = str(new_read_cov)
        f.write('\t'.join(line) + '\n')
    
print("Down-sampled to {} reads.".format(str(total_out_read)))
