#!/bin/bash

# this code removes the line breaks in fasta files

infasta=$1
outfasta=$2

awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' $infasta > $outfasta