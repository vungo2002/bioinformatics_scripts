#!/usr/bin/env python3
from sys import argv
import os 

vcf_file = argv[1]
outfile = argv[2]

cmd1 = """bcftools query -i'FILTER="PASS"' -i'GT="1/1" || GT="1|1"' -f'%CHROM  %POS     %ID      %REF     %ALT     %QUAL    %FILTER  %INFO    %FORMAT\n' {}| wc -l > {}""".format(vcf_file, outfile)
print(cmd1)
os.system(cmd1)
cmd2 = """bcftools query -i'FILTER="PASS"' -i'GT="0/1" || GT="0|1" || GT="1|0"' -f'%CHROM  %POS     %ID      %REF     %ALT     %QUAL    %FILTER  %INFO    %FORMAT\n' {}| wc -l >> {}""".format(vcf_file, outfile)
print(cmd2)
os.system(cmd2)

lines = open(outfile).read().strip().split()
hom_count = int(lines[0])
het_count = int(lines[1])

with open(outfile,'w') as f:
    f.write("hetcount,hom_count,ratio\n")
    line = str(het_count) + ',' + str(hom_count) + ',' + str(het_count/hom_count)
    f.write(line)
    print(line)
