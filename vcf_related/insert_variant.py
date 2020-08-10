from sys import argv

"""
Usage: 
infile = argv[1] # phased VCF file called by Beagle (separated by '|' instead of '/')
reffile = argv[2]  # genomic reference, fasta file (eg hg19, hg38, etc)
chrom = argv[3] # chromosome to insert, only works per chromosome, not whole genome
target_sample = argv[4] # target sample to look at in the VCF file
hap_idx = int(argv[5]) # usually 0 or 1
outfile = argv[6] # an altered genome in FASTA format
"""


def read_phased_vcf_data(infile):
    """ input: the phased vcf file produced by Beagle 
        output: an list of coordinates and genotypes from the vcf file
    """
    
    data = []
    for line in open(infile):
        if line.startswith("#CHROM"):
            header = line.strip().split()
            #print (header)

        if line.startswith("#") == False:
            #print (line)
            tmp = line.strip().split()
            data += [tmp]
    return header,data

def get_haplotypes(data, target_sample):   
    target_idx = -1
    for i in range(len(header)):
        if header[i] == target_sample:
            target_idx = i

    haplotypes = []
    for i in range(len(data)):
        chrom = data[i][0]
        loc = data[i][1]
        ref = data[i][3]
        alts = data[i][4]
        gt = data[i][target_idx]
        if gt != "0|0": 
            haplotypes += [[chrom, loc, ref, alts, gt]]
    return haplotypes

def load_ref(fastafile, target_chrom=None):
    """ load the reference genome from a single fasta file"""
    chroms = open(fastafile).read().split(">")[1:]
    genome = {}
    for chrom in chroms:
        tmp = chrom.split()
        chr_name = tmp[0]
        if target_chrom != None:
            if chr_name != target_chrom:
                continue
        seq = ''.join(tmp[1:])
        genome[chr_name] = seq
    return genome

def insert_variants(ref_seq, haplotypes, hap_idx):
    """ input: choose to insert either haplotype 0 or 1 to the reference
        output: the altered reference
    """
    # copy the seqs first, also convert the string into list of characters
    new_seq = ''
    
    print ("total numver of variants", len(haplotypes))
    counter = 0
    
    prev_var_idx = 0 # tracks the location of the previous variant (VCF coordinate is 1-based)
    total_added = 0 # tracks the difference in bp number after altering 
    for line in haplotypes:
        if counter%10000 == 0:
            print(counter)
        counter += 1
        

        chrom = line[0]
        loc_s = int(line[1]) - 1
        ref = line[2]
        loc_e = loc_s + len(ref)
        alts = [ref] + line[3].split(',')
        alt_idx = int(line[4].split("|")[hap_idx])
        alt = alts[alt_idx].replace("*",'')


        new_seq += ref_seq[prev_var_idx:loc_s]
        #print(ref, alt)
        new_seq = new_seq + alt
        prev_var_idx = loc_e
        
        total_added += len(alt) - len(ref)
    new_seq += ref_seq [prev_var_idx:]
        
    print("total bases dif", total_added)
    return new_seq



### Main program

infile = argv[1] # phased VCF file called by Beagle
ref_file = argv[2]  # genomic reference, fasta file
chrom = argv[3] # chromosome to insert
target_sample = argv[4] # target sample to look at in the VCF file
hap_idx = int(argv[5])
outfile = argv[6] 


header, data = read_phased_vcf_data(infile)
haplotypes = get_haplotypes(data, target_sample)
del data

#load genome 
genome = load_ref(ref_file, chrom)
# alter the chrom
altered_chrom = insert_variants(genome[chrom], haplotypes, hap_idx)

out = open(outfile, 'w')
out.write(">" + chrom + '_altered\n' )
out.write(altered_chrom)
out.close()