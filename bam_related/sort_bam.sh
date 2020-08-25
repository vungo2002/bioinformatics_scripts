#!/bin/bash 

inbam=$1
cpus=$2
outbam=$3

samtools sort -@ $cpus $inbam > $outbam 