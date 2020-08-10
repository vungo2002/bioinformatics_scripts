# this script takes 2 bam files and mix them together at the ratio specified 
# this assumes that two bam files are in equal coverages


from sys import argv
import os

bam1 = argv[1] # fetal/placental bam file
ff1 = float(argv[2]) # the fraction of the bam1 that will be in the output
bam2 = argv[3] # bam file 2, maternal healthy
ff2 = 1.0 - ff1 # fraction of bam2 that will be in the output
out = argv[4]


# get a fraction of bam1
seed = 1111
fraction = str(seed + ff1)
subsampled_bam1 = bam1 +'.subsample.'+str(ff1)+'.bam'
cmd = "samtools view -s %s %s -b> %s" % (fraction, bam1, subsampled_bam1)
print(cmd)
os.system(cmd)


# get a fraction of bam2
#fraction = str(seed + ff2)
#subsampled_bam2 = bam2 +'.subsample.'+str(ff2)+'.bam'
#cmd = "samtools view -s %s %s -b> %s" % (fraction, bam2, subsampled_bam2)
#print(cmd)
#os.system(cmd)

# merging
#cmd = "samtools merge -f -@ 2 %s %s %s" %(out, subsampled_bam1, subsampled_bam2)
cmd = "samtools merge -f -@ 2 %s %s %s" %(out, subsampled_bam1, bam2)
print(cmd)
os.system(cmd)

# sorting the merged file 
srtbam = out+'.srt.tmp'
cmd = "samtools sort -@ 2 -m 2500M %s > %s" %(out, srtbam)
print(cmd)
os.system(cmd)

# remove tmp file
cmd = "rm %s" %(out)
print(cmd)

# change the name of srtbam to out
cmd = "mv %s %s" %(srtbam, out)
print(cmd)
os.system(cmd)

# index the result
cmd = "samtools index %s" %(out)
print(cmd)
os.system(cmd)


# cleaning tmp files
cmd = "rm %s" %(subsampled_bam1)
print(cmd)
os.system(cmd)
#cmd = "rm %s" %(subsampled_bam2)
#print(cmd)
#os.system(cmd)

print("Done.")


