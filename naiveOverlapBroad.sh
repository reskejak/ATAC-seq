#!/bin/bash

##########################

# naiveOverlapBroad.sh
# Jake Reske
# Michigan State University, 2019
# reskejak@msu.edu
# https://github.com/reskejak

# Computing ENCODE-defined naive overlapping broadPeak set from two ATAC or ChIP individual replicate broadPeak and replicate-pooled broadPeak sets
# Definition from Kundaje et al.: "Find pooled peaks that overlap Rep1 and Rep2 where overlap is defined as the fractional overlap wrt any one of the overlapping peak pairs >= 0.5".

# usage: bash naiveOverlapBroad.sh treat1_peaks.broadPeak treat2_peaks.broadPeak treat_pool_peaks.broadPeak > treat_overlap_peaks.broadPeak
# ensure to firstly make executable by: chmod u+x naiveOverlapBroad.sh

# dependencies: bedtools (for intersectBed / bedtools intersect)

# can extend function to accept more replicates as desired

##########################

REP1=$1
REP2=$2
POOL=$3

# broadPeak
intersectBed -wo \
-a ${POOL} -b ${REP1} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | \
cut -f 1-9 | sort | uniq | \
intersectBed -wo \
-a stdin -b ${REP2} | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | \
cut -f 1-9 | sort | uniq
