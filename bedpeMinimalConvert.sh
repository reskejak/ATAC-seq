#!/bin/bash

##########################

# bedpeMinimalConvert.sh
# Jake Reske
# Michigan State University, 2019
# reskejak@msu.edu
# https://github.com/reskejak

# convert standard bedtools BEDPE format to "minimal" BEDPE format accepted by MACS2 for peak calling
# note, can take 20-40 minutes or longer per large (>10 million reads) BEDPE file

# usage: bash bedpeMinimalConvert.sh treat1.bedpe > treat1.minimal.bedpe
# ensure to firstly make executable by: chmod u+x bedpeMinimalConvert.sh

##########################

FILE=$1

# iterate through columns 2, 3, 5 and 6. Update the min and the max values and print in "minimal" BEDPE format
awk '{min=max=$2 ; for (i=3 ; i <= 6; i++) if (i!=4) {($i<min)?min=$i:($i>max)?max=$i:min}};{printf $1"\t%s\t%s\t"$7"\n", min, max}' $FILE -
