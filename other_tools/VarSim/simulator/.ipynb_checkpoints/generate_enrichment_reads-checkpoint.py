#!/usr/bin/env python3

from sys import argv
import os
import shutil


def load_chrom_seq(genome_fasta, chrom):
    """
    This function only loads 1 sequence from the genome fasta file
    """
    seq = ''
    inseq = False
    for line in open(genome_fasta):
        if line.strip() == '>' + chrom:
            # this is inside the seq
            inseq = True
            continue
        elif line.startswith(">"): # this is a different chromosome
            inseq = False
        if inseq:
            seq += line.strip()
    genome = {}
    genome[chrom] = seq
    return genome 

def insert_variant_2(ref_seq, locs, refs_to_be_replaced, alts_to_replace): 
    """
    Input: a sequence, a list of loci on that sequence, and a list of things to convert the ref into 
    Output: altered sequence, old_locs, new_locs 

    - old_locs are the old locs that was in the input
    - new_locs are the old coordinates since the alternation can change the length of the sequence thus shift many locations
    """
    new_seq = ''
    counter = 0

    new_coords = [] # this stores the coordinates of the inserted variants on the newly altered chromosome
    old_coords = [] # this stores old coordinates of inserted variants

    prev_var_coord = 0 # tracks the location of the end of the previous variant (VCF coordinate is 1-based), this is on the reference genome
    total_added = 0 # tracks the difference in bp number after altering chrom

    for i in range(len(locs)): 
        loc = locs[i]
        counter += 1
        if counter % 10000 == 0:
            print(counter)

        loc_s = loc # this function did not need to convert 1-based to 0-based because it assumes 0-based
        ref = refs_to_be_replaced[i]
        loc_e = loc_s + len(ref) # loc_s and loc_e are coords on the reference genome
        alt = alts_to_replace[i]

        new_seq += ref_seq[prev_var_coord:loc_s] # adding the space between the current variant and the previous one
        new_coords += [len(new_seq)] 
        old_coords += [loc_s]
        
        new_seq = new_seq + alt
        prev_var_coord = loc_e
        
        total_added += len(alt) - len(ref)
    new_seq += ref_seq [prev_var_coord:]
        
    print("inserted number of variants", counter)
    print("total bases dif", total_added)
    
    return new_seq, old_coords, new_coords

def make_art_sim_reads(fastafile, fold_cov, output_prefix, read_length=100, fragment_size_mean=200, fragment_size_std=20, art_illumina_PATH='/art_sim/'):
    cmd = art_illumina_PATH + "art_illumina -ss HS25 -na -i %s -p -f %s -l %s -m %s -s %s -o %s" %(fastafile, fold_cov, read_length, fragment_size_mean, fragment_size_std, output_prefix)
    #print(cmd)
    os.system(cmd)
    return




# TODO: add more ART related metrics


### MAIN program
inputfile = argv[1] # "example_input.bed"
ref_fasta = argv[2] # the genomic fasta file (hg38, hg19, etc.) "/efs/data/genomes/hg38/hg38.fa" 

inputs = {}
for line in open(inputfile):
    tmp = line.strip().split()
    #print(tmp)
    loc_name = tmp[0] +'_' + tmp[1] + '-' + tmp[2]
    inputs[loc_name] = tmp 
inputs.keys()


alt_out_fasta_dir = "alt_seqs/"
ref_out_fasta_dir = "ref_seqs/"
try:
    os.mkdir(alt_out_fasta_dir)
except FileExistsError:
    shutil.rmtree(alt_out_fasta_dir)
    os.mkdir(alt_out_fasta_dir)
try:
    os.mkdir(ref_out_fasta_dir)
except FileExistsError:
    shutil.rmtree(ref_out_fasta_dir)
    os.mkdir(ref_out_fasta_dir)


# Generate the fastas sequences for the regions of interest
chromosomes = ['chr' + str(x) for x in range(1,23)] +['chrX', 'chrY']
print("Making the fasta sequences...")
for chrom in chromosomes:
    locs = []
    refs_to_be_replaced = []
    alts_to_replace = []
    regions_of_interest = []
    targets = []
    for loc in inputs:
        line = inputs[loc]
        if line[0] != chrom:
            continue
        else:
            locs += [int(line[1])]
            refs_to_be_replaced += [line[4]]
            alts_to_replace += [line[5]]
            regions_of_interest += [[int(line[9]), int(line[10])]]
            targets += [loc]
    if len(locs) == 0: # nothing to do
        print("no locs to alter in", chrom)
        continue
    else:
        print("loading...", chrom)
        ref_seq = load_chrom_seq(ref_fasta, chrom)[chrom]
        new_seq, old_coords, new_coords = insert_variant_2(ref_seq, locs, refs_to_be_replaced, alts_to_replace)
        # recalculate the regions of interest based on the shifts between new_coords and old_coords
        for i in range(len(old_coords)):
            shift = new_coords[i] - old_coords[i]
            new_region = [regions_of_interest[i][0] + shift, regions_of_interest[i][1] + shift]

            # now extract the sequences out and put them into a fasta file
            seq_name = targets[i]
            alt_out = open(alt_out_fasta_dir + seq_name + '.alt.fa', 'w')
            alt_out.write('>' + targets[i] +'-alt' +'\n')
            alt_out.write(new_seq[new_region[0]:new_region[1]] +'\n')
            alt_out.close()
            
            ref_out = open(ref_out_fasta_dir + seq_name +'.ref.fa', 'w')
            ref_out.write('>' + targets[i] +'-ref' +'\n')
            ref_out.write(ref_seq[regions_of_interest[i][0]:regions_of_interest[i][1]] +'\n')
            ref_out.close()


# Generate the reads with ART

# do the alt_seqs first
for file in os.listdir(alt_out_fasta_dir):
    if file.endswith('.fa') == False:
        continue
    #print(file)
    loc_name = file.split('.')[0]
    line = inputs[loc_name]
    #print(line)
    fold_cov = float(line[6])*float(line[7])
    fasta_file = alt_out_fasta_dir + file
    output_prefix = fasta_file.replace('.fa','.art.')
    make_art_sim_reads(fasta_file, fold_cov, output_prefix)

# do the ref_seqs second
for file in os.listdir(ref_out_fasta_dir):
    if file.endswith('.fa') == False:
        continue
    #print(file)
    loc_name = file.split('.')[0]
    line = inputs[loc_name]
    fold_cov = (1-float(line[6]))*float(line[7])
    fasta_file = ref_out_fasta_dir + file
    output_prefix = fasta_file.replace('.fa','.art.')
    make_art_sim_reads(fasta_file, fold_cov, output_prefix)



# in the alt_seqs/ 
cmds = []
cmds += ["cat %s > %s" %(alt_out_fasta_dir+'/*.1.fq', 'alt_reads.1.fq')]
cmds += ["cat %s > %s" %(alt_out_fasta_dir+'/*.2.fq', 'alt_reads.2.fq')]
# in the ref_seqs/
cmds += ["cat %s > %s" %(ref_out_fasta_dir+'/*.1.fq', 'ref_reads.1.fq')]
cmds += ["cat %s > %s" %(ref_out_fasta_dir+'/*.2.fq', 'ref_reads.2.fq')]
cmds += ['cat alt_reads.1.fq ref_reads.1.fq > enrichment_combined.sim.reads.1.fq']
cmds += ['cat alt_reads.2.fq ref_reads.2.fq > enrichment_combined.sim.reads.2.fq']

for cmd in cmds:
    print(cmd)
    os.system(cmd)
# clean up
files_to_remove = ['alt_reads.1.fq', 'alt_reads.2.fq','ref_reads.1.fq', 'ref_reads.2.fq']
for file in files_to_remove:
    os.remove(file)

dirs_to_remove = ['alt_seqs/', 'ref_seqs/']
for dir in dirs_to_remove:
    shutil.rmtree(dir)