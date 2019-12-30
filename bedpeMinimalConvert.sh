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

# define coordSort function
function coordSort()
{
while read line; do
      ary=(${line})
      for((i=0; i!=${#ary[@]}; ++i)); do
         for((j=i+1; j!=${#ary[@]}; ++j)); do
              if (( ary[i] > ary[j] )); then
                    ((key=ary[i]))
                    ((ary[i]=ary[j]))
                    ((ary[j]=key))
              fi
         done
      done
      echo ${ary[@]}
done < ${1}
}

# get paired-end coordinates from BEDPE file
awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\n",$2,$3,$5,$6}' $FILE > coords.$FILE

# sort coordinates and arrange in minimal format
coordSort coords.$FILE | paste - $FILE | awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\n",$5,$1,$4,$11}' -

# remove intermediate coordinates file
rm coords.$FILE
