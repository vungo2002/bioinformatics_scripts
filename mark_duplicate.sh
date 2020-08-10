#!/bin/bash

inbam=$1
outbam=$2
metric_file=$3

gatk MarkDuplicates -I $inbam -O $outbam -M $metric_file
