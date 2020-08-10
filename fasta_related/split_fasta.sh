#!/bin/bash

fafile=$1 # absolute path
outdir=$2

mkdir $outdir
cd $outdir

csplit -s -z $fafile '/>/' '{*}'
for i in xx* ; do \
    n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
    mv "$i" "$n.fa" ; \
    done