# this script parses the output directory of A DRAGEN pipeline and returns certain metrics in json format

import os 
from sys import argv
import json 

indir = argv[1]
outfile = argv[2]


keys = ["Uniformity of coverage (PCT > 0.2*mean) over target region", 
        "Uniformity of coverage (PCT > 0.2*mean) over genome",
        "Het/Hom ratio",
        "Aligned reads,",
        "Aligned reads in target region"]



json_keys_dict = {"Uniformity of coverage (PCT > 0.2*mean) over target region":'uniformity_of_coverage_pct_gt_02mean_over_target_region',
                    "Uniformity of coverage (PCT > 0.2*mean) over genome":"uniformity_of_coverage_pct_gt_02mean_over_genome",
                    "Het/Hom ratio":"variants_het_to_hom_ratio_pass",
                    "Aligned reads,":"aligned_reads",
                    "Aligned reads in target region":"aligned_reads_in_target_region"
                    }


results = {}
for file in os.listdir(indir):
    if file.endswith(".csv"):
        #print(file)
        for line in open(indir + '/' + file):
            for k in keys:
                if k in line:
                    tmp = line.strip().split(',')
                    if k != "Aligned reads in target region":
                        #print(json_keys_dict[k],tmp)
                        results[json_keys_dict[k]] = float(tmp[-1])
                    else:
                        #print(k, line)
                        results[json_keys_dict[k]] = float(tmp[-2])
                        results["aligned_reads_in_target_region_pct"] = float(tmp[-1])

                    
print(results)
with open(outfile,'w') as f:
    f.write(json.dumps(results, indent=2))
    