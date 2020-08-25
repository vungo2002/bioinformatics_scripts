#!/bin/bash

bam1=$1
bam2=$2
outbam=$3
cpus=$4
samtools merge -r -@ $cpus $outbam $bam1 $bam2 