# this script takes the bedtools coverage -d (single bp) input and output the average coverages 
# per chromosome  # specify these later


from sys import argv
import os


infile = argv[1] # bed cov file
print ("input is", infile)

tmp_file = infile.replace(".cov",".tmp")
cmd = "cut -f1,5 %s | sort|uniq -c > %s" %(infile, tmp_file)
print("COMMAND is", cmd)
os.system(cmd)

# spit the result into chromosomes

results = {}
for line in open(tmp_file):
    tmp = line.split()
    chrom = tmp[1]
    if chrom not in results:
        results[chrom] = [tmp]
    else:
        results[chrom] += [tmp]

for chrom in results:
    total_bp = 0
    read_covs = 0
    for line in results[chrom]:
        total_bp += int(line[0])
        read_covs += int(line[0])*int(line[2])
    avg_cov = float(read_covs)/total_bp
    print(chrom, avg_cov)




