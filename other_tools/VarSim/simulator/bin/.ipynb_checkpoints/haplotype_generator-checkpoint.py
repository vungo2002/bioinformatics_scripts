#!/usr/bin/env python3

import subprocess
import os
import argparse
from argparse import RawTextHelpFormatter
import logging 
import sys

# TODO: add function to assess input file type when loading hap_gen_input
# Make sure that the BED files are SORTED

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


def get_roi_seqs(bedFile, refFasta):
    """
    1. From region of interest BED file, output FASTA file containig
    sequences for each region.
    2. Load sequences for each region
    """
    tmpFAFile = bedFile + '.tmpfasta'
    subprocess.run(['bedtools', 'getfasta', '-fi', refFasta, '-bed', bedFile, '-fo', tmpFAFile])
    refSeqs = load_ref(tmpFAFile)
    os.remove(tmpFAFile)

    return refSeqs


def convert_roi_names_to_coordinates(listOfROINames):
    listOfROIs = []
    for roi in listOfROINames:
        tmp = roi.split(":")
        tmp = [tmp[0]] + [int(x) for x in tmp[1].split('-')]
        listOfROIs += [tmp]
    return listOfROIs


def assign_to_roi(loc, listOfROIs, one_based=False):
    """
    From a location, find a region of interest that contain said location and output new loc.
    INPUT:
        - location in the form of [CHROMOSOME, POS] # note, POS is 0-based
        - ROIs in the form of [[CHROMOSOME, START, END]] # START and END are also 0-based
    OUTPUT:
        - ROI that contains the input location
        - New cordinates with respect to that ROI in the form of [ROI, START, END]
    """
    for roi in listOfROIs:
        if roi[0] == loc[0]:
            if one_based:
                loc[1] = loc[1] - 1
            if roi[1] <= loc[1] and loc[1] < roi[2]:
                newPos = loc[1] - roi[1]
                roiName = roi[0] + ':' + str(roi[1]) + '-' + str(roi[2])
                return roiName, newPos
    return None


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


def convert_hapGenInput_to_new_coords(hapGenInput, listOfROIs, one_based=False):
    """
    Convert the hapGenInput to new coordinates based on the reference sequences

    If the locations of variants are 1-based, use one_based = True.
    """
    newHapGenInput = []
    for line in hapGenInput:
        loc = line[:2]
        try:
            roi, newPos = assign_to_roi(loc, listOfROIs, one_based=one_based)
        except TypeError:  # the return is None, the variant is not inside any ROI
            continue
        newHapGenInput += [[roi, newPos] + line[2:]]

    return newHapGenInput


def insert_variants(refSeqs, seqName, hapGenInput, one_based=False):
    # ref_seq is the seqeuence of the chrom of interest
    # seq_name is the name of the ref_seq
    # haplotype_calls is in the form ['chr1', '30001230', 'G', 'A']

    # copy the seqs first, also convert the string into list of characters
    newSeq = ''
    counter = 0
    new_coords = []  # this stores the coordinates of the inserted variants on the newly altered chromosome
    old_coords = []  # this stores old coordinates of inserted variants
    prev_var_coord = 0  # tracks the location of the end of the previous variant (VCF coordinate is 1-based), this is on the reference genome
    total_added = 0  # tracks the difference in bp number after altering chrom
    for line in hapGenInput:  # $haplotypes is in the form [[chrom, loc, ref, alts, gt]]
        # print(line)
        chrom = line[0]
        if chrom != seqName:
            continue
        counter += 1
        if counter % 10000 == 0:
            print(counter)

        if one_based:
            loc_s = int(line[1]) - 1  # convert 1-based to 0-based
        else:
            loc_s = int(line[1])

        ref = line[2]
        loc_e = loc_s + len(ref)  # loc_s and loc_e are coords on the reference genome
        alt = line[3].replace("-", '').replace("*", '')  # if '-', it's a deletion
        newSeq += refSeqs[prev_var_coord:loc_s]  # adding the space between the current variant and the previous one
        new_coords += [len(newSeq)]
        old_coords += [loc_s]
        newSeq = newSeq + alt
        prev_var_coord = loc_e
        total_added += len(alt) - len(ref)
    newSeq += refSeqs[prev_var_coord:]

    # print("Inserted number of variants", counter)
    # print("Total bases dif", total_added)

    return newSeq, old_coords, new_coords

def print_converted_input(newHapGenInput, conInFile, tag=''):
    """
    Print out converted hap_gen_input file 
    """
    with open(conInFile, 'w') as f:
        for i in range(len(newHapGenInput)):
            tmp = newHapGenInput[i]
            if tag != '':
                tmp[0] = tmp[0] + tag
            line = '\t'.join([str(x) for x in tmp])
            f.write(line +'\n')
    return

def generate_haplotype_seqs(hapGenInput, refSeqs, one_based=False, tag=''):
    """
    INPUT: Haplotype Input, reference sequences
    OUTPUT: Haplotype sequences in FASTA format
    """
    haploSeqs = {}
    for seq in refSeqs:
        newSeq, old_coords, new_coords = insert_variants(refSeqs[seq], seq, hapGenInput, one_based=one_based)
        haploSeqs[seq + tag] = newSeq
        #print(seq, new_coords)
    return haploSeqs


def main():
    """
    This script is just to create haplotype based on the requirements
    """
    global VERSION
    VERSION = "1.0.0"
    
    args = get_arguments()
    
    hapGenBed = args.input
    refFa = args.ref
    outFa = args.outFasta
    bedROI = args.bedROI
    conInFile = args.conInFile
    isOneBased = args.isOneBased
    tag = args.tag
    
    if outFa == 'simHaplotype.fa':
        logging.info("output not named, using default")
    if refFa == None:
        logging.error("Missing reference. Exit.")
        sys.exit(2)
    if hapGenBed == None:
        logging.error("Missing input. Exit.")
        sys.exit(1)
    if bedROI == None:
        logging.error("Missing regions-of-interest BED file. Exit.")
        sys.exit(3)
    if tag != None:
        tag = '_' + tag
    else:
        tag = ''
        
    hapGenInput = load_hap_gen_input(hapGenBed)
    refSeqs = get_roi_seqs(bedROI, refFa)
    listOfROIs = convert_roi_names_to_coordinates(list(refSeqs.keys()))
    newHapGenInput = convert_hapGenInput_to_new_coords(hapGenInput, listOfROIs, one_based=isOneBased)
    haploSeqs = generate_haplotype_seqs(newHapGenInput, refSeqs, one_based=isOneBased, tag=tag)
    
    with open(outFa, 'w') as f:
        for seq in haploSeqs:
            f.write('>' + seq + '\n')
            f.write(haploSeqs[seq] + '\n')
            
    print_converted_input(newHapGenInput, conInFile, tag)

    return



def get_arguments():


    parser = argparse.ArgumentParser(prog='haplotype_generator.py',
                                    description="""
    Author: Vu Ngo
    This program takes:
    - haplotype_input_bed file
    - reference fasta sequence
    - regions of interest
    Output:
    - fasta sequences for each region of interest with variants inserted""", formatter_class=RawTextHelpFormatter)

    parser.add_argument('-i','--input', action='store', dest='input',
                       help="REQUIRED. Input containing locations and variants to insert to reference.")
    parser.add_argument('-r','--ref', action='store', dest='ref',
                       help="REQUIRED. Reference FASTA file that matches with the coordinates in the input.")
    parser.add_argument('-b', '--roi-bed', action='store', dest='bedROI',
                       help="REQUIRED. BED file containing regions of interest to be extracted from the simulated genome haplotype.")
    parser.add_argument('-o','--out-fasta', action='store', dest='outFasta', default="simHaplotype.fa",
                        help='The output haplotype FASTA name, if not specified, a default file name "simHaplotype.fa" will be used')
    parser.add_argument('-c','--converted', action='store', dest='conInFile', default='convertedInput.txt',
                        help='Converted input file that has coordinates match with the new FASTA file.')
    parser.add_argument('--one-based', action='store', dest="isOneBased", default=False,
                       help="BOOLEAN: TRUE or FALSE. This specifies whether the input is in 0-based or 1-based coordinate. Default = FALSE")
    parser.add_argument('-t','--tag', action='store', dest='tag', default=None,
                       help="A tag to add to each of the simulated haplotype sequences' names. Default is None")
    
    parser.add_argument('--version', '-v', action='version', version="%(prog)s "  + VERSION)

    args = parser.parse_args()
    return args

if __name__=="__main__":
    main()