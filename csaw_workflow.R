# csaw_workflow.R

# csaw workflow for ATAC-seq differential accessibility analysis, in R
# Jake Reske
# Michigan State University, 2020
# reskejak@msu.edu
# https://github.com/reskejak

# Machine-readable version of Figure 6 workflow for ATAC-seq DA analysis with csaw
# will describe methods for using 1) pre-defined peaks from MACS2 as well as 2) csaw de novo enriched window calling by local enrichment, 
# and normalization methods including 1) TMM on binned counts and 2) loess-based for trended biases

# brief aside: use gc() to help clear memory after intensive commands if crashes/errors occur

# example experimental design: n=2 mouse ATAC-seq biological replicates for two conditions: treat and control 

library(GenomicRanges)
library(csaw)

########################################
########################################
########################################

# starting from MACS2 filtered broadpeaks

# read replicate broadPeak files
treat1.peaks <- read.table("treat1_peaks.filt.broadPeak", sep="\t")[,1:3]
treat2.peaks <- read.table("treat2_peaks.filt.broadPeak", sep="\t")[,1:3]
control1.peaks <- read.table("control1_peaks.filt.broadPeak", sep="\t")[,1:3]
control2.peaks <- read.table("control2_peaks.filt.broadPeak", sep="\t")[,1:3]
colnames(treat1.peaks) <- c("chrom", "start", "end")
colnames(treat2.peaks) <- c("chrom", "start", "end")
colnames(control1.peaks) <- c("chrom", "start", "end")
colnames(control2.peaks) <- c("chrom", "start", "end")

# read naive overlap broadPeak files
treat.overlap.peaks <- read.table("treat_overlap_peaks.filt.broadPeak", sep="\t")[,1:3]
control.overlap.peaks <- read.table("control_overlap_peaks.filt.broadPeak", sep="\t")[,1:3]
colnames(treat.overlap.peaks) <- c("chrom", "start", "end")
colnames(control.overlap.peaks) <- c("chrom", "start", "end")

# convert to GRanges objects
treat1.peaks <- GRanges(treat1.peaks)
treat2.peaks <- GRanges(treat2.peaks)
treatn.peaks <- GRanges(treatn.peaks)
treat.overlap.peaks <- GRanges(treat.overlap.peaks)
control1.peaks <- GRanges(control1.peaks)
control2.peaks <- GRanges(control2.peaks)
controln.peaks <- GRanges(controln.peaks)
control.overlap.peaks <- GRanges(control.overlap.peaks)

# define consensus peakset

# one method: union of all replicate peak sets for both conditions
treat.peaks <- union(treat1.peaks, treat2.peaks)
control.peaks <- union(control1.peaks, control2.peaks)
all.peaks <- union(treat.peaks, control.peaks)

# another method: intersect between biological replicates; union between both experimental conditions
treat.peaks <- intersect(treat1.peaks, treat2.peaks)
control.peaks <- intersect(control1.peaks, control2.peaks)
all.peaks <- union(treat.peaks, control.peaks)

# yet another method: union between naive overlapping peak sets
all.peaks <- union(treat.overlap.peaks, control.overlap.peaks)

##############################
# specify paired-end BAMs
pe.bams <- c("control1.sorted.noDups.filt.noMT.bam", "control2.sorted.noDups.filt.noMT.bam",
	     "treat1.sorted.noDups.filt.noMT.bam", "treat2.sorted.noDups.filt.noMT.bam")

##############################
# read mm10 blacklist
blacklist <- read.table("~/ref_genome/mm10.blacklist.bed", sep="\t")
colnames(blacklist) <- c("chrom", "start", "end")
blacklist <- GRanges(blacklist)

# define read parameters
standard.chr <- paste0("chr", c(1:19, "X", "Y")) # only use standard chromosomes
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)

##############################
# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)

##############################
# MACS2 peaks only: filter low abundance peaks
library("edgeR")
peak.abundances <- aveLogCPM(asDGEList(peak.counts)) 
peak.counts.filt <- peak.counts[peak.abundances > -3, ] # only use peaks logCPM > -3
# few or no peaks should be removed; modify as desired

##############################

# get paired-end fragment size distribution
control1.pe.sizes <- getPESizes(control1.pe.bam)
control2.pe.sizes <- getPESizes(control2.pe.bam)
treat1.pe.sizes <- getPESizes(treat1.pe.bam)
treat2.pe.sizes <- getPESizes(treat2.pe.bam)
gc()
# plot
hist(treat1.pe.sizes$sizes) # repeat for all replicates and conditions

# for analysis with csaw de novo enriched query windows, select a window size that is greater than the majority of fragments

##############################
# count BAM reads in, e.g. 300 bp windows
counts <- windowCounts(pe.bams, width=300, param=param) # set width as desired from the fragment length distribution analyses

# filter uninteresting features (windows) by local enrichment
# local background estimator: 2kb neighborhood
neighbor <- suppressWarnings(resize(rowRanges(counts), width=2000, fix="center")) # change width parameter as desired
wider <- regionCounts(pe.bams, regions=neighbor, param=param) # count reads in neighborhoods
filter.stat <- filterWindows(counts, wider, type="local") 
counts.local.filt <- counts[filter.stat$filter > log2(3),] # threshold of 3-fold increase in enrichment over 2kb neighborhood abundance; change as desired

###############################
# count BAM background bins (for TMM normalization)
binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)

##########################################
# NORMALIZATION

# method 1: MACS2 peaks only, TMM normalization based on binned counts
peak.counts.tmm <- peak.counts.filt
peak.counts.tmm <- normFactors(binned, se.out=peak.counts.tmm)

# method 2: MACS2 peaks only, csaw loess-normalization
peak.counts.loess <- peak.counts.filt
peak.counts.loess <- normOffsets(peak.counts.loess, type="loess", se.out=TRUE)
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

# method 3: csaw de novo peaks by local enrichment, TMM normalization based on binned counts
counts.local.tmm <- counts.local.filt
counts.local.tmm <- normFactors(binned, se.out=counts.local.tmm)

# method 4: csaw de novo peaks by local enrichment, csaw loess-normalization
counts.local.loess <- counts.local.filt
counts.local.loess <- normOffsets(counts.local.loess, type="loess", se.out=TRUE)
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

#########################################
# DIFFERENTIAL ACCESSIBILITY ANALYSIS

# set working windows for the desired analysis
working.windows <- peak.counts.tmm # MACS2 peaks only, standard TMM normalization based on binned counts
# working.windows <- peak.counts.loess # MACS2 peaks only, for trended biases
# working.windows <- counts.local.tmm # csaw de novo peaks by local enrichment, standard TMM normalization based on binned counts
# working.windows <- counts.local.loess # csaw de novo peaks by local enrichment, for trended biases
# SEE THE CSAW MANUAL FOR MORE INFO ON NORMALIZATION METHODS
###########

# setup design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)
colnames(y$counts) <- c("control1", "control2", "treat1", "treat2")
rownames(y$samples) <- c("control1", "control2", "treat1", "treat2")
y$samples$group <- c("control", "control", "treat", "treat")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("treat", "control")
# design

# stabilize dispersion estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
results <- glmQLFTest(fit, contrast=makeContrasts(treat-control, levels=design))
# head(results$table)

# combine GRanges rowdata with DA statistics
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# should merge some peaks; change as desired

# use most significant window as statistical representation for p-value and FDR for merged windows
tab.best <- getBestTest(merged.peaks$id, results$table)
head(tab.best)
# combine merged peaks window range with statistics
final.merged.peaks <- merged.peaks$region
final.merged.peaks@elementMetadata <- cbind(final.merged.peaks@elementMetadata, tab.best[,-1])
final.merged.peaks <- final.merged.peaks[order(final.merged.peaks@elementMetadata$FDR), ] # sort by FDR
final.merged.peaks # all windows

# filter by FDR threshold
FDR.thresh <- 0.05 # set as desired
final.merged.peaks.sig <- final.merged.peaks[final.merged.peaks@elementMetadata$FDR < FDR.thresh, ]
final.merged.peaks.sig # significant differentially-accessible windows

write.table(final.merged.peaks, "treat_vs_control_csaw_DA-windows_all.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(final.merged.peaks.sig, "treat_vs_control_csaw_DA-windows_significant.txt", sep="\t", quote=F, col.names=T, row.names=F)

###########################################

# Generate MA plot
library(ggplot2)

final.merged.peaks$sig <- "n.s."
final.merged.peaks$sig[final.merged.peaks$FDR < FDR.thresh] <- "significant"

ggplot(data=data.frame(final.merged.peaks),
       aes(x = logCPM, y = logFC, col = factor(sig, levels=c("n.s.", "significant")))) + 
  geom_point() + scale_color_manual(values = c("black", "red")) + 
  geom_smooth(inherit.aes=F, aes(x = logCPM, y = logFC), method = "loess") + # smoothed loess fit; can add span=0.5 to reduce computation load/time
  geom_hline(yintercept = 0) + labs(col = NULL)
