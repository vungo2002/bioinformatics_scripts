#!/bin/bash
inputpdf=$1
output_prefix=$2
pdftoppm $inputpdf $output_prefix -png
