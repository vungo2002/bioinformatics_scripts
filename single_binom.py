from scipy import stats
import numpy as np
from sys import argv


def reduce_num(vector, targetsum):
    """ input a vector of ints, reduce the number of those so that the elements sum up to targetsum,
    """
    cursum = sum(vector)
    ratio = targetsum/cursum
    
    return [int(x*ratio) for x in vector]

def get_inform_loci_binom(read_profile, pval_thres = 0.0001):
    """
    This function identifies informative loci by calculating the binomial p-val of observing a read distribution 
    """
    het_p = 0.5
    ref_p = 0.999
    alt_p = 1 - ref_p

    inform_loci = []
    targetsum = 100

    for index in range(len(read_profile)):
        i = read_profile[index]
        if len(i) > 2:  # ignore cases where there are more than 2 alleles for a locus
            continue
        if sum(i) > targetsum:   # reduce the number of reads to reduce library prep bias's effect
            i = reduce_num(i, targetsum)
        try: 
            ratio = i[0]/sum(i)
        except ZeroDivisionError:
            continue
        if ratio < 0.1:
            pval = stats.binom_test(x = i, p = alt_p)
        elif ratio > 0.9:
            pval = stats.binom_test(x = i, p = ref_p)
        else:
            pval = stats.binom_test(x = i, p = het_p)

        if pval < pval_thres:
            #print(i, largest_pval, labels[largest_index])
            inform_loci += [index]

    return inform_loci


def calc_minorF(read_profile, informative_loci):
    minor_Fs = []
    for i in informative_loci:
        rc = read_profile[i]
        factor_1 = min(rc)*2
        factor_2 = sum(rc) - factor_1
        total = sum(rc)
        try:
            minor_f = min(factor_1, factor_2)/total
        except ZeroDivisionError:
            continue

        minor_Fs += [minor_f]
    if len(minor_Fs) == 0:
        return "Na"
    return np.median(minor_Fs)


### MAIN Program

input_file = argv[1]


depth_thres = 2
# read the read profile 
read_profile = []
for line in open(input_file):
    if line.startswith("#"):
        continue
    tmp = line.strip().split()
    ad = tmp[-1].split(":")[1].split(",")
    ad = [int(x) for x in ad]
    for num in ad:
        if num < depth_thres:
            continue
    #print(ad)
    read_profile += [ad]

informative_loci = get_inform_loci_binom(read_profile, pval_thres = 0.0001)

ff_estimated = calc_minorF(read_profile, informative_loci)

print ("Number of informative loci: " + str(len(informative_loci)))
print("FF: "+ str(ff_estimated))