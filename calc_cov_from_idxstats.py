from sys import argv
import os

infile = argv[1] # bam file that has .bai file with it 

tmpfile = infile.replace("bam","tmp")
cmd = "samtools idxstats %s | grep chr|grep -v '_' > %s" %(infile, tmpfile)
os.system(cmd)

for line in open(tmpfile):
    tmp = line.strip().split()
    print('\t'.join(tmp) + '\t' +str(int(tmp[2])*150/float(tmp[1])))

os.system("rm "+ tmpfile)