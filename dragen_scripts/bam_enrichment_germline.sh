#!/bin/bash 

inbam=$1
targetbed=$2
outdir=$3

base=$(basename $inbam)
prefix=${base/.bam}


mkdir $outdir
/opt/edico/bin/dragen -f -r /usr/local/illumina/genomes/hg38_alt_aware/DRAGEN/8 -b $inbam --enable-variant-caller true --output-directory $outdir --output-file-prefix  $prefix --vc-emit-ref-confidence GVCF --vc-target-bed $targetbed  --vc-target-bed-padding 150 