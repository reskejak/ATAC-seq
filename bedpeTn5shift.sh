#!/bin/bash

##########################

# bedpeTn5shift.sh
# Jake Reske
# Michigan State University, 2019
# reskejak@msu.edu
# https://github.com/reskejak

# intended to shift BEDPE coordinates +4 and -5 bp, respectively, to compensate for Tn5 insertion in ATAC data as described by Buenrostro et al. 2013 (PMID: 24097267).

# usage: bash bedpeTn5shift.sh treat1.bedpe > treat1.tn5.bedpe
# ensure to firstly make executable by: chmod u+x bedpeTn5shift.sh

##########################

FILE=$1

awk -F $'\t' 'BEGIN {OFS = FS}{ \
if ($9 == "+") {$2 = $2 + 4; $6 = $6 - 5} \
else if ($9 == "-") {$3 = $3 - 5; $5 = $5 + 4} \
print $0}' $FILE
