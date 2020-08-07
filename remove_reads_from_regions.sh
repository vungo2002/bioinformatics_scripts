inbam=$1
bed_region_to_remove=$2
outbam=$3


bedtools intersect -abam $inbam -b $bed_region_to_remove -v > $outbam
