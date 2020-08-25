# this script takes a bam file, look at the different samples in the header, then change all of them to one sample (specified in the input)
import os 
from sys import argv
import sys

inbam = argv[1]
new_sample_name = argv[2]
outbam = argv[3]

if " " in new_sample_name or "\t" in new_sample_name:
    print("Invalid sample name")
    sys.exit(0)


# get the header of the bam file into a tmp file 
header_file = inbam + ".header"
cmd = "samtools view -H %s > %s " %(inbam, header_file) 
os.system(cmd)

# read the header and replace the sample name with the new sample name 

new_header_file = open(inbam +".new_header", "w")
for line in open(header_file):
    if line.startswith("@RG"):
        tmp = line.split('\t')
        for i in range(len(tmp)):
            ele = tmp[i]
            if ele.startswith("SM"):
                tmp[i] = "SM:" + new_sample_name
        newline = '\t'.join(tmp)
        new_header_file.write(newline)
    else:
        new_header_file.write(line)
new_header_file.close()

# now fix the original bam file with a new header 

cmd = "samtools reheader %s %s > %s" %(inbam +".new_header", inbam, outbam)
os.system(cmd)

# remove tmp files 
os.system("rm " + header_file)
os.system("rm " + inbam +".new_header")
 