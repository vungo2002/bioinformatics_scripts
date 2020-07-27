

"""This program will take 2 BAM files, a coverage N, and and normalize BAM1 down to coverage N using the read distribution of BAM2

USAGE:
python3 normBam.py bam1 bam2 covN bam1_normed chrom_exceptions
"""

from sys import argv
import sys
import pysam
import os
import time
from collections import defaultdict
import random as rd



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


# read the results from the idxstats file 
def read_idxstats(idxstatsfile):
    """ read the idxstats results (grep chr|grep -v '_' """
    covs = {}  # each entry contain the cov, the length of the region, the number of reads in the region, respectively
    readlength = 150
    total_read_BP = 0
    total_length = 0
    for line in open(idxstatsfile):
        tmp = line.strip().split()
        #print (tmp)
        chrom = tmp[0]
        readbp = int(tmp[2])*readlength
        region_length = int(tmp[1])
        cov = readbp/float(region_length)
        covs[chrom] = [cov, int(region_length), int(tmp[2])]

        total_read_BP += readbp
        total_length += region_length

    avg_cov = total_read_BP/float(total_length)
    return covs, avg_cov




### MAIN PROGRAM ###


bam1 = argv[1] #"../data/August_Illumina_WGS/1X_BAM_wo_dups/20-P9_S1.0.592.ds.nodup.ds0.052.1X.bam"
bam2 = argv[2]
target_cov = float(argv[3])
bam1_normed = argv[4]
chrom_exception_file = argv[5] # a file listing chromosomes that will not be normalized 


# reading BAM 1

indexfile = bam1+'.bai'
print(indexfile)
# check to see if the index file is there
if os.path.exists( indexfile):
    print ("exists")
    
else:
    print("error: no index file found")
    sys.exit() 
    
# get the idxstats results from the index file 
idxstatsfile = indexfile.replace(".bai",".idxstats")
cmd = "samtools idxstats %s |grep chr | grep -v '_' > %s" % (bam1, idxstatsfile)
print(cmd)
os.system(cmd)

covs, avg_cov = read_idxstats(idxstatsfile)
print(covs)

# reading BAM 2


indexfile = bam2+'.bai'
print(indexfile)
# check to see if the index file is there
if os.path.exists( indexfile):
    print ("exists")
    
else:
    print("error: no index file found")
    sys.exit() 
    
# get the idxstats results from the index file 
idxstatsfile = indexfile.replace(".bai",".idxstats")
cmd = "samtools idxstats %s |grep chr | grep -v '_' > %s" % (bam2, idxstatsfile)
print(cmd)
os.system(cmd)

covs2, avg_cov2 = read_idxstats(idxstatsfile)


# #determine the number of reads needed for a certain average coverage

# use covs2 as referrence:
print (avg_cov2)
ds_factor_bam2 = target_cov/avg_cov2

target_num_reads = {} # this stores the number of reads per chromosome required for BAM2 to get down to desired cov

for chrom in covs2:
    target_num_reads[chrom] = ds_factor_bam2*covs2[chrom][2]
    #print(chrom,covs2[chrom], ds_factor_bam2*covs2[chrom][2])
    
target_ds_factors = {} # this stores the ds factor to downsize BAM1 down to desired cov
for chrom in target_num_reads:
    target_ds_factors[chrom] = target_num_reads[chrom]/(covs[chrom][2] + 0.001)
    

# listing chroms that will not be normalized relatively to chromosome 1
no_norm_list = [] 
for line in open(chrom_exception_file):
    no_norm_list += [line.strip()]

for chrom in no_norm_list:
    print ("no norming chrom", chrom)
    target_ds_factors[chrom] = target_ds_factors['chr1']
# do not norm ChrY to the same as ChrX 
target_ds_factors['chrY'] = target_ds_factors['chr1'] # this will keep the relative ratio the same no matter what
target_ds_factors['chrX'] = target_ds_factors['chr1']
target_ds_factors['chrEBV'] = target_ds_factors['chr1']
target_ds_factors['chrM'] = target_ds_factors['chr1']

print("chrom", "cov_before", "ds_factor", "cov after norm will be")
for chrom in target_num_reads:
    print(chrom, covs[chrom][0],target_ds_factors[chrom], covs[chrom][0]*target_ds_factors[chrom])


for chrom in target_ds_factors: # input check to see if valid
    if target_ds_factors[chrom] > 1.0:
        print ("ERROR: fraction of BAM 1 is larger than 1.0, please change the inputs to make it work")
        sys.exit()



print ("Subsampling BAM1 to normalize the read distribution... ")

start_time = time.time()

tmpBAM1 = bam1 +".tmp"
samfile = pysam.AlignmentFile(bam1, "rb")
pairedreads = pysam.AlignmentFile(tmpBAM1, "wb", template=samfile)

for chrom in target_ds_factors:
    counter = 0
    print(chrom, "ds factor",target_ds_factors[chrom])
    for read1, read2 in read_pair_generator(samfile,chrom):
        if rd.random() <= target_ds_factors[chrom]:
            counter += 1
            pairedreads.write(read1)
            pairedreads.write(read2)
    print( counter)

pairedreads.close()
samfile.close()

print("--- %s seconds ---" % (time.time() - start_time))   


print ("sorting the result")
start_time = time.time()
cmd = "samtools sort %s > %s" % (tmpBAM1, bam1_normed)
os.system(cmd)
print("--- %s seconds ---" % (time.time() - start_time))   

# remove the tmp file
os.system("rm %s" %(tmpBAM1))

print ("indexing result file")
cmd = "samtools index %s" %(bam1_normed)
os.system(cmd)

print("Done.")


