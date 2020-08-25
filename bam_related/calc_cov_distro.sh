#!/bin/bash

inbam=$1
samtools index $inbam
bamCoverage -b $inbam --outFileFormat bedgraph --ignoreDuplicates --binSize 50000 -o ${inbam/.bam/.cov.bedgraph}
