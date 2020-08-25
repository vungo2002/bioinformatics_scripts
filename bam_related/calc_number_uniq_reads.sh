inbam=$1
outcount=$2

samtools view -c -F 4 -F 1024 -F 2048 $inbam > $outcount
