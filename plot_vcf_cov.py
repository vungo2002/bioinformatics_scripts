#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from sys import argv
import os 

def readGZVCF(infile, ad_index):
    """
    This function reads a vcf.gz file and reports the location-genotype-AD from the file 
    """
    # extract the file using vcftools then grep just the loc + the genotype info
    # puts it in a temporary file

    tmpfile = infile+'.tmp'
    if os.path.exists(tmpfile) == False:

        cmd1 = "bcftools view %s |grep -v '#'|grep -v 'chrM'|cut -f1-5,10 > %s " %(infile, tmpfile)
        #print(cmd1)
        os.system(cmd1)
    
    # read the temp file and load the info into a dictionary 

    dict = {}
    for line in open(tmpfile):
        tmp = line.strip().split('\t')
        loc = ':'.join(tmp[:2])
        alt_allele = tmp[4].split(',')
        if len(alt_allele)>1:
            continue
        else:
            ADs = [int(x) for x in tmp[-1].split(":")[ad_index].split(',')]
            dict[loc] = ADs 

    return dict



## MAIN
infile = argv[1]
AD_index = int(argv[2])
vcf_dict = readGZVCF(infile, AD_index)
cov_thres = 1

covs = []
for loc in vcf_dict:
    cov = sum(vcf_dict[loc])
    if cov < cov_thres:
        continue
    else:
        covs += [cov]

plt.hist(covs,bins=1000)
#plt.xlim(0,100)
plt.title(infile.split("/")[-1] + " coverages" )
plt.savefig(infile.replace(".vcf.gz",'.cov.png'))
plt.close()


covs_lognormed= [np.log(x) for x in covs] 
plt.hist(covs_lognormed, bins = 100)
plt.xlim(0,6)
plt.title(infile.split("/")[-1] + " log-normalized coverages")
plt.savefig(infile.replace(".vcf.gz",'.logcov.png'))
plt.close()
