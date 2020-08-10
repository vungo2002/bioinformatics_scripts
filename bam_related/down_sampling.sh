inbam=$1
ds_factor=$2
outbam=$3

gatk DownsampleSam -I $inbam -O $outbam -P $ds_factor -R null