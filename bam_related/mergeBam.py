import os
from sys import argv

"""
Input: 2 bam files
Output: 1 output bam file

"""

bam1 = argv[1]
bam2 = argv[2]
outbam = argv[3]


# merge the bam files 
print("merging")
cmd = "samtools merge -@ 2 %s %s %s" %(outbam, bam1, bam2)
print(cmd)
os.system(cmd)

# sort the new bam file
print("sorting...")
srtbam = outbam.replace(".bam",".srt.bam")
cmd = "samtools sort -@ 2 -m 3000M %s > %s" %(outbam, srtbam)
print(cmd)
os.system(cmd)
# index the new bam file
print("indexing...")
cmd = "samtools index %s" %(srtbam)
print(cmd)
os.system(cmd)
print("done")