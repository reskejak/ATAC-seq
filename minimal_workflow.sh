#!/bin/bash

# ATAC-seq workflow and software commands
# March 2018
# Michigan State University
# reskejak@msu.edu

DIR=("/mnt/home/usr/foo")
REFDIR=("/mnt/home/usr/ref_genome") # location of reference genomes
BIN=("/mnt/home/usr/bin")
samples=("foo_treat1 foo_treat2 foo_treatn foo_control1 foo_control2 foo_controln")

#########################

cd ${DIR}

# remove excess filename architecture, paired end reads
for i in $samples;
do mv ${i}_L000_R1_001.fastq.gz ${i}_R1.fastq.gz;
mv ${i}_L000_R2_001.fastq.gz ${i}_R2.fastq.gz;
done

#########################
# fastqc on raw files
for i in $samples;
do fastqc ${i}_R1.fastq.gz;
fastqc ${i}_R2.fastq.gz;
done

##########################
# make trimmed output directory
mkdir trimmed

# trim galore
for i in $samples;
do trim_galore --paired -o trimmed/ ${i}_R1.fastq.gz ${i}_R2.fastq.gz;
done

cd trimmed

##########################
# fastqc on trimmed files
for i in $samples;
do fastqc ${i}_R1_trimmed.fq.gz;
fastqc ${i}_R2_trimmed.fq.gz
done

###########################
# bowtie2 for alignment

module load bowtie2
module load samtools

# build mm10 index
cd ${REFDIR}
bowtie2-build mm10.fa mm10_bt2_index

# align
# note: alignment only takes ~20 minutes with 64 gb RAM and 16 procs
# supply two hours job time to be safe for unzipping, aligning, zipping, and sorting bam
cd ${DIR}/trimmed
for i in $samples;
do gunzip ${i}_R1_trimmed.fq.gz;
gunzip ${i}_R2_trimmed.fq.gz;
bowtie2 -p 16 --very-sensitive \
-X 1000 \
-x ${REFDIR}/mm10_bt2_index \
-1 ${i}_R1_trimmed.fq \
-2 ${i}_R2_trimmed.fq | samtools view -bS - > ${i}.bam;
gzip ${i}_R1_trimmed.fq;
gzip ${i}_R2_trimmed.fq;
done
#note: -X 1000 indicates no alignments over 1000 bp will be collected

#samtools sort bam
samtools sort -T /tmp/${i}.sorted -o ${i}.sorted.bam ${i}.bam

###############
# remove mitochondrial reads

HARVARD=("${BIN}/ATAC-seq") # location of Harvard ATAC-seq module
module load samtools

for i in $samples;
do samtools view -h ${i}.sorted.bam | python ${HARVARD}/atacseq/removeChrom.py - - chrM | samtools view -b - > ${i}.noMT.sorted.bam;
done

#sort by coordinates again
samtools sort -T /tmp/${i}.sorted.noMT -o ${i}.sorted.noMT.bam ${i}.noMT.sorted.bam
# wc -l after sorting at this point yields identical result as removeChrom output report

########################
# remove PCR duplicates

PICARD=("${BIN}/picard") # location of picard tools

for i in $samples;
do java -jar ${PICARD}/build/libs/picard.jar MarkDuplicates \
I=${i}.sorted.noMT.bam \
O=${i}.noDups.sorted.noMT.bam \
M=${i}_dups.txt \
REMOVE_DUPLICATES=true;
done
# command takes ~15-20 minutes/sample in interactive dev node

######################
# sort by read name and update SAM flags following alignment post-processing
# and convert to bedpe
# https://www.biostars.org/p/149119/

module load samtools
module load bedtools

# loop for sort by read name, fix mate flags, and convert properly-paired reads to bedpe
for i in $samples;
do samtools sort -n ${i}.noDups.sorted.noMT.bam ${i}.sorted.noDups.noMT;
samtools fixmate ${i}.sorted.noDups.noMT.bam ${i}.fixed.bam;
samtools view -bf 0x2 ${i}.fixed.bam | bedtools bamtobed -i stdin -bedpe > ${i}.fixed.bedpe;
done

# example for calculating mapping statistics
# samtools flagstat ${i}.fixed.bam

# note: bedpe corresponds to exactly half of "properly paired" due to read collapse

########################
# TN5 transposase shifting of BEDPE files

for i in $samples;
do awk -F $'\t' 'BEGIN {OFS = FS}{ \
if ($9 == "+") {$2 = $2 + 4; $6 = $6 - 5} \
else if ($9 == "-") {$3 = $3 - 5; $5 = $5 + 4} \
print $0}' ${i}.fixed.bedpe > ${i}.tn5.bedpe;
done

#######################
# convert files to properly formatted "minimal" bedpe for macs2
# macs2 simplified bedpe format is: [chr] [start] [stop] [name]

for i in $samples;
do awk 'BEGIN{OFS="\t"}{printf "%s\t%s\t%s\t%s\n",$2,$3,$5,$6}' ${i}.tn5.bedpe > ${i}.coords.bed;
bash sort.sh ${i}.coords.bed | paste - ${i}.tn5.bedpe | awk 'BEGIN{OFS="\t"}\
{printf "%s\t%s\t%s\t%s\n",$5,$1,$4,$11}' - > ${i}.final.tn5.bedpe;
rm ${i}.coords.bed;
done
# takes ~20-40 minutes per sample

####################
# sort.sh
# sort each row horizontally
# usage: bash sort.sh coords.bed > coords.sorted.bed
#####################
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
####################
# https://www.unix.com/shell-programming-and-scripting/180835-sort-each-row-horizontally-awk-any.html

#########################
# pool biological replicates

POOLED=("foo_treat foo_control") # group identifier, i.e. experimental condition

for i in $POOLED;
do cat ${i}1.final.tn5.bedpe ${i}2.final.tn5.bedpe > ${i}.cat.final.tn5.bedpe;
done

# short computation time <2 minutes each sample

################################
# Peak calling by macs2
# NOTE: Ben Johnson suggests calling broadpeaks for ATAC-seq

# export bin to PATH for calling macs2
export PATH=/mnt/home/usr/.local/bin:$PATH

ALLSAMPLES=("foo_treat1 foo_treat2 foo_treatn foo_treat.cat foo_control1 foo_control2 foo_controln foo_control.cat")
# real samples and pooled replicates

# narrow peaks
for i in $ALLSAMPLES;
do macs2 callpeak -t ${i}.final.tn5.bedpe \
-f BEDPE \
-n ${i}_narrow \
-g mm \
--keep-dup all;
done
# broad peaks
for i in $ALLSAMPLES;
do macs2 callpeak -t ${i}.final.tn5.bedpe \
-f BEDPE \
-n ${i}_broad \
-g mm \
--broad \
--broad-cutoff 0.05 \
--keep-dup all;
done

###################################
# mask repeats by blacklist filtering
module load bedtools

for i in $ALL_SAMPLES;
do bedtools intersect -v -a ${i}_peaks.narrowPeak \
-b mm10.blacklist.bed | grep -P 'chr[\dXY]+[ \t]' | awk 'BEGIN{OFS="\t"} {print $0}' > ${i}_peaks.filt.narrowPeak;
done

##################################
# Naive overlap thresholding for MACS2 peak calls
# (not required for DiffBind differential-peak calling)
module load bedtools

for i in $POOLED;
# Find pooled broadPeaks that overlap Rep1 and Rep2
# overlap is defined as the fractionaloverlap wrt any one of the overlapping peak pairs  >= 0.5
do intersectBed -wo \
-a ${i}.cat_broad_peaks.broadPeak -b ${i}1_broad_peaks.broadPeak | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq | \
intersectBed -wo \
-a stdin -b ${i}2_broad_peaks.broadPeak | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$12-$11; if (($19/s1 >= 0.5) || ($19/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq > ${i}.PooledInRep1AndRep2.broadPeak;
done

for i in $POOLED;
# Find pooled narrowPeaks that overlap Rep1 and Rep2
# overlap is defined as the fractionaloverlap wrt any one of the overlapping peak pairs  >= 0.5
do intersectBed -wo \
-a ${i}.cat_narrow_peaks.narrowPeak -b ${i}1_narrow_peaks.narrowPeak | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq | \
intersectBed -wo \
-a stdin -b ${i}2_narrow_peaks.narrowPeak | \
awk 'BEGIN{FS="\t";OFS="\t"}{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {print $0}}' | \
cut -f 1-10 | sort | uniq > ${i}.PooledInRep1AndRep2.narrowPeak;
done

#############################
# Differential analysis between 2 conditions
# https://groups.google.com/forum/#!topic/idr-discuss/Q2B2biJ9is0

# DiffBind in R, via HPCC

# load latest version of R
module swap GNU GNU/4.9
module load OpenMPI/1.10.0
module load R/3.3.2
R

# set local R library
.libPaths("/mnt/home/usr/R/local/R_libs")

setwd("/mnt/home/usr/foo")

library("DiffBind")

# import sample file csv
# see DiffBind manual for generation of this file
samples <- read.csv("foo-sample-file.csv")

# construct experiment DBA object, DESeq2 analysis
foo<- dba(minOverlap = 2, 
	sampleSheet = "foo-sample-file.csv", 
	peakCaller = "bed", 
	config=data.frame(AnalysisMethod=DBA_DESEQ2, fragmentSize=150))
# for fragment sizes, use macs2 .xls output file param
# for example...
foo$config$fragmentSize <- c(129, 140, 137, 134, 102, 91, 88)

# count reads in open-chromatin intervals
foo <- dba.count(foo)

# setup contrasts for differential-chromatin analysis
foo <- dba.contrast(foo, minMembers=2, categories=DBA_FACTOR)
	
# differential-chromatin analysis
foo <- dba.analyze(foo, bFullLibrarySize=FALSE)
# DESeq2 only works for my data when bFullLibrarySize=FALSE (default is =TRUE)
# Look into this DiffBind parameter for more information about usage...

# summarize output
foo.db <- dba.report(foo)

# write output to csv
write.csv(foo.db, file="foo_DESeq2-diff-peaks.csv")

# plot sample correlation by peaks
pdf(plot(foo, contrast=1), file="foo-samples-peaks.pdf")

#################################
# HOMER - peak annotation and motif finding

HOMER=("${BIN}/homer/bin") # location of HOMER tools
export PATH=$PATH:${BIN}/homer/bin

# peak annotation, gene ontology, and genome ontology
annotatePeaks.pl foo_DESeq2-diff-peaks_new.bed mm10 \
-size given \
-annStats KOvsWT_DESeq2_Homer-annotation-stats.txt \
-go KOvsWT_DESeq2_HOMER-GO \
-genomeOntology KOvsWT_DESeq2_HOMER-genomeOntology > KOvsWT_DESeq2_HOMER-annotation.txt

# peak annotation
annotatePeaks.pl foo_DESeq2-diff-peaks.txt mm10 -size given > HOMER_anno_size-given.txt # must input TSV without header column
# output does not change with or without "-size given" flag

# find motifs
findMotifs.pl foo_DESeq2-diff-peaks.txt mouse foo_MotifResults -start -400 -end 100 -len 8,10 -p 4 
# "This will search for motifs of length 8 and 10 from -400 to +100 relative to the TSS, using 4 threads (i.e. 4 CPUs)"
# be wary of 

#################################
# creating UCSC browser tracks
# based on Tao Liu's bdgcmp method
# https://github.com/taoliu/MACS/wiki/Build-Signal-Track#Run_MACS2_bdgcmp_to_generate_foldenrichment_and_logLR_track

# re-run macs2 (broad peaks) with bedGraph flag for each true replicate
# use -B --SPMR for normalizing to fragments per million reads
for i in $samples;
do macs2 callpeak -t ${i}.final.tn5.bedpe \
-f BEDPE \
-n ${i}_broad_norm \
-g mm \
-B \
--SPMR \
--broad \
--broad-cutoff 0.05 \
--keep-dup all;
done

# use macs2 bdgcmp to build track
# description of bdgcmp options for generating signal to noise statistic visualization
# https://groups.google.com/forum/#!topic/macs-announcement/yefHwueKbiY
# log likelihood ratio is very clean; could also try fold enrichment (see documentation) but exhibits much more noise
# see documentation for -p flag (pseudocount) addition to avoid log(0) calculation
for i in $samples;
do macs2 bdgcmp -t ${i}_broad_norm_treat_pileup.bdg \
-c ${i}_broad_norm_control_lambda.bdg \
-o ${i}_broad_norm_logLR.bdg \
-m logLR \
-p 0.00001;
done

module load bedtools
# bdg2bw, from Tao Liu (MACS2)
# https://gist.github.com/taoliu/2469050
for i in $samples;
do ${BIN}/bdg2bw ${i}_broad_norm_logLR.bdg ${REFDIR}/mm10.chrom.sizes;
done

# wigCorrelate biological replicates for statistic
for i in $POOLED;
do ${BIN}/wigCorrelate ${POOLED}1_broad_norm_logLR.bdg ${POOLED}2_broad_norm_logLR.bdg;
done

module load bedtools
# convert bed to bigWig in one pipe command
# http://seqanswers.com/forums/showthread.php?t=4731
for i in $samples;
do bedtools sort -i ${i}_broad_peaks.filt.broadPeak | \
bedtools genomecov -bg -i stdin -g ${REFDIR}/mm10.chrom.sizes | \
wigToBigWig -clip stdin ${REFDIR}/mm10.chrom.sizes ${i}_broad_peaks.filt.bigWig;
done

# then upload and generate public URL for bigWig files and create UCSC track hub for visualization
# for more information: https://genome.ucsc.edu/goldenpath/help/hgTrackHubHelp.html
