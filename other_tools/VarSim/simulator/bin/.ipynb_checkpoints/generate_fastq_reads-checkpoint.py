#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import logging
import argparse
from argparse import RawTextHelpFormatter

logging.basicConfig(level=logging.DEBUG)

def load_ref(fastaFile, targetChrom=None):
    """
    Load the reference genome from a single fasta file.
    - If sequence name is specified, only load that sequence.
    """
    chroms = open(fastaFile).read().split(">")[1:]
    refSeqs = {}
    for chrom in chroms:
        tmp = chrom.split()
        chr_name = tmp[0]
        if targetChrom != None:
            if chr_name != targetChrom:
                continue
        seq = ''.join(tmp[1:])
        refSeqs[chr_name] = seq
    return refSeqs

def load_hap_gen_input(hapGenBed):
    """
    Read haplotype generator input BED file
    """
    hapGenInput = []
    for line in open(hapGenBed, 'r'):
        line = line.strip().split('\t')
        if len(line) > 3:
            line[1] = int(line[1])
            hapGenInput += [line]
    return hapGenInput



def make_art_sim_reads(ART_PATH, read1_prof, read2_prof, fasta_file, fold_cov, output_prefix, **kwargs):
    
    read_length=kwargs['read_length']
    frag_size_mean=kwargs['frag_size_mean']
    frag_size_std=kwargs['frag_size_std']
    
    cmd = [ART_PATH, '-1', read1_prof, '-2', read2_prof, '-na', '-i', fasta_file, '-p', '-f', fold_cov, '-l', read_length, 
           '-m', frag_size_mean, '-s', frag_size_std, '-o', output_prefix]
    #logging.info(' '.join(cmd))
    try: 
        subprocess.run(cmd, check=True, capture_output=True)
    except:
        logging.error(' '.join(cmd))
    return


def make_art_sim_reads_all_regions(ART_PATH, hap_gen_input, hap_seqs, read1_prof, read2_prof, fastq_prefix, tmp_out_dir, **kwargs):
    """
    INPUT: hap_gen_input, hap_seqs, read1_profile, read2_profile

    Generates fastq files for each region, the combine them all together into 2 FASTQ files

    NOTES:
    - If there are multiple variants in a region/sequence, the coverage for that region will be the average of those variants.
    """

    read_length=kwargs['read_length']
    frag_size_mean=kwargs['frag_size_mean']
    frag_size_std=kwargs['frag_size_std']
    
    try:
        os.makedirs(tmp_out_dir)
    except FileExistsError:
        logging.warning("Dir exists:" + tmp_out_dir)
    for seq in hap_seqs:
        # look for all variants in this seq, there isn't any, continue
        covs = []
        for ele in hap_gen_input:
            if ele[0] == seq:
                covs += [float(ele[4])]
        if covs == []: #this sequence doesn't contain any variant
            #print(seq + ' does not contain any variant')
            continue
        else:
            # make a fasta file for just this sequence
            fasta_file = 'tmp.' + seq + '.fa'
            with open(fasta_file, 'w') as f:
                f.write('>' + seq + '\n')
                f.write(hap_seqs[seq] + '\n')
            fold_cov = str(np.mean(covs))
            output_prefix = tmp_out_dir + '/' + seq + "_"
            make_art_sim_reads(ART_PATH, read1_prof, read2_prof, fasta_file, fold_cov, output_prefix,
                              read_length=read_length, frag_size_mean=frag_size_mean, frag_size_std=frag_size_std)
            os.remove(fasta_file) # remove tmp file
    
    # combine all output .fq files by read number
    subprocess.run(["cat " + tmp_out_dir + "/*_1.fq > " + fastq_prefix + "_R1.fq"], shell=True, check=True, capture_output=True)
    subprocess.run(["cat " + tmp_out_dir + "/*_2.fq > " + fastq_prefix + "_R2.fq"], shell=True, check=True, capture_output=True)
    
    # remove the individual fastq files
    subprocess.run(["rm -r " + tmp_out_dir], shell=True, check=True, capture_output=True)
    return

def main():
    global VERSION
    VERSION = "1.0.0"

    args = get_arguments()
    ART_PATH = args.ART_PATH
    hap_gen_input_file = args.input
    hap_seqs_file = args.hap_seqs
    read1_prof = args.read1_prof
    read2_prof = args.read2_prof
    out_prefix = args.out_prefix
    read_length = args.read_length
    frag_size_mean = args.frag_size_mean
    frag_size_std = args.frag_size_std


    hap_gen_input = load_hap_gen_input(hap_gen_input_file)
    hap_seqs = load_ref(hap_seqs_file)
    tmp_out_dir = out_prefix +'_tmp'

    make_art_sim_reads_all_regions(ART_PATH, hap_gen_input, hap_seqs, read1_prof, read2_prof, out_prefix, tmp_out_dir,
                                   read_length=read_length, frag_size_mean=frag_size_mean, frag_size_std=frag_size_std)
    return

def get_arguments():
    parser = argparse.ArgumentParser(prog='generate_fastq_reads.py',
                                    description="""
    Author: Vu Ngo
    This program takes:
    - haplotype_input_bed file
    - haplotype fasta file
    - read 1 and read 2 ART error model
    - ART_PATH
    - output_prefix
    Output:
    - fastq reads in the form of output_prefix_R[1,2].fq""", formatter_class=RawTextHelpFormatter)
    parser.add_argument('-a', '--art-path', action='store', dest='ART_PATH',
                        help="REQUIRED. Executabe path for art_illumina")
    parser.add_argument('-1', '--read1', action='store', dest='read1_prof',
                        help="REQUIRED. Error model for read1.")
    parser.add_argument('-2', '--read2', action='store', dest='read2_prof',
                        help="REQUIRED. Error model for read2.")
    parser.add_argument('-i','--input', action='store', dest='input',
                        help="REQUIRED. haplotype_input_bed.")
    parser.add_argument('-r','--ref', action='store', dest='hap_seqs',
                        help="REQUIRED. Haplotype FASTA file that matches with the coordinates in the input.")
    parser.add_argument('-o','--out-prefix', action='store', dest='out_prefix',
                        help='REQUIRED. Output prefix for FASTQ files.')
    parser.add_argument("-l","--read_length", action='store', dest='read_length', default='100',
                        help="length of read. Default is 100")
    parser.add_argument('-m', '--frag-size-mean', action='store', dest='frag_size_mean', default="200",
                       help="mean fragment size. Default is 200.")
    parser.add_argument('-s', '--frag-size-std', action='store', dest='frag_size_std', default='20',
                       help="standard deviation of fragment size. Default is 20.")
    parser.add_argument('--version', '-v', action='version', version="%(prog)s "  + VERSION)

    args = parser.parse_args()
    return args

if __name__=="__main__":
    main()