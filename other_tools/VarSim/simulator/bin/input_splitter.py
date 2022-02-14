#!/usr/bin/env python3

import os
import subprocess
import logging
import sys
import numpy as np
import argparse
from argparse import RawTextHelpFormatter

"""
This script is to split the main input into separate haplotype inputs. 
The results will be fed into haplotype generator
"""

VERSION = 'v1.0.0'
logging.basicConfig(level=logging.DEBUG)

def load_hap_fraction(hap_frac_str):
    tmp = hap_frac_str.strip().split(":")[1].split('\t')
    labels = tmp[0].split(';')
    values = tmp[1].split(';')
    fracs = {}
    for i in range(len(labels)):
        fracs[labels[i]] = float(values[i])
    return fracs

def header_dict(header_list):
    """
    Return a dict with the index of each header in the list"""
    header_idx = {}
    for i in range(len(header_list)):
        header_idx[header_list[i]] = i
    return header_idx

def main_input_split(main_input_file, out_dir, noise=0.0):
    """
    - Takes main input file, output directory, noise level (optional)
    - Split into hap_gen_inputs accordingly
    - For each hap, gives a number of read coverage (on average, they follow the HAPLOTYPE_FRACTION metric)
    """
    # create output dir
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass
    # read the fractions of haplotypes
    for line in open(main_input_file, 'r'):
        if line.startswith("#HAPLOTYPE_FRACTION"):
            fracs = load_hap_fraction(line)

    # check to see if hap_input.txt files already exist, if so, remove them
    hap_input_files = []
    for hap in fracs:
        f = out_dir +'/' + hap + '_input.txt'
        if os.path.exists(f):
            logging.info("Haplotype input file {} already exist. It will be replaced.".format(f))
            os.remove(f)
        hap_input_files.append(f)

    # read the REF + ALT alleles and assign them to each HAPLOTYPE
    for line in open(main_input_file, 'r'):
        if line.startswith('#CHROM'):
            header = line.strip().replace("#",'').split('\t')
            header_idx = header_dict(header)
        elif line.startswith("#") == False:
            tmp = line.strip().split('\t')
            total_cov = int(tmp[header_idx['COV']])
            ref = tmp[header_idx["REF"]]
            alt = tmp[header_idx["ALT"]]
            alleles = [ref] + alt.split(',')
            #print(total_cov, ref, alt, alleles)
            for hap in fracs:
                fout = out_dir +'/' + hap+'_input.txt'
                with open(fout, 'a') as f:
                    line = tmp[:3]
                    line += [alleles[int(tmp[header_idx[hap]])]]
                    hap_cov = str(total_cov*fracs[hap] * max(0, np.random.normal(loc=1.0, scale=noise))) # add some noise to hap_cov
                    line += [hap_cov]
                    f.write('\t'.join(line)+'\n')
                    #print(line)

    return fracs, hap_input_files


def main():
    args = get_arguments()
    mainInput = args.input
    out_dir = args.out_dir
    noise = float(args.noise)
    main_input_split(mainInput, out_dir, noise=noise)
    return

def get_arguments():

    parser = argparse.ArgumentParser(prog='haplotype_generator.py',
                                    description="""
    Author: Vu Ngo
    This program takes:
    - mnv_simulator main input file
    - output directory

    Output:
    - haplotype_generator input files
    - desired coverages for each region """, formatter_class=RawTextHelpFormatter)

    parser.add_argument('-i','--input', action='store', dest='input',
                       help="REQUIRED. Main MNV_simulator input file.")
    parser.add_argument('-o','--out-dir', action='store', dest='out_dir',
                       help="REQUIRED. Output directory. Haplotype_generator input files will be stored here.")
    parser.add_argument('-n','--noise', action='store', dest='noise', default=0.0,
                       help="Noise level, between 0 and 1.0. Default is 0.0.")
    parser.add_argument('--version', '-v', action='version', version="%(prog)s "  + VERSION)
    args = parser.parse_args()
    return args

if __name__=="__main__":
    main()