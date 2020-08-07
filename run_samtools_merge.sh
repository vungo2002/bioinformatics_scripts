inbam_1=$1
inbam_2=$2
outbam=$3

samtools merge -r $outbam $inbam_1 $inbam_2