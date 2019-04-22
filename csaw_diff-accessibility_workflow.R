# in R
# csaw workflow for ATAC differential-accessibility
# April 2019

# will describe methods for using 1) pre-defined peaks from MACS2 as well as 2) de novo csaw peak calling by local enrichment, 
# and normalization methods including 1) TMM on binned counts and 2) loess-based for trended biases

# brief aside: use gc() to help clear memory after intensive commands if crashes/errors occur

########################################
# starting from MACS2 filtered broadpeaks

library("GenomicRanges")

# read broadPeak files
treat1.peaks <- read.table("treat1_broad_peaks.filt.broadPeak", sep="\t")[,1:3]
treat2.peaks <- read.table("treat2_broad_peaks.filt.broadPeak", sep="\t")[,1:3]
treatn.peaks <- read.table("treatn_broad_peaks.filt.broadPeak", sep="\t")[,1:3]
control1.peaks <- read.table("control1_broad_peaks.filt.broadPeak", sep="\t")[,1:3]
control2.peaks <- read.table("control2_broad_peaks.filt.broadPeak", sep="\t")[,1:3]
controln.peaks <- read.table("controln_broad_peaks.filt.broadPeak", sep="\t")[,1:3]
colnames(treat1.peaks) <- c("chrom", "start", "end")
colnames(treat2.peaks) <- c("chrom", "start", "end")
colnames(treatn.peaks) <- c("chrom", "start", "end")
colnames(control1.peaks) <- c("chrom", "start", "end")
colnames(control2.peaks) <- c("chrom", "start", "end")
colnames(controln.peaks) <- c("chrom", "start", "end")

# convert to GRanges objects
treat1.peaks <- GRanges(treat1.peaks)
treat2.peaks <- GRanges(treat2.peaks)
treatn.peaks <- GRanges(treatn.peaks)
control1.peaks <- GRanges(control1.peaks)
control2.peaks <- GRanges(control2.peaks)
controln.peaks <- GRanges(controln.peaks)

# define consensus peakset
# intersect between biological replicates; union between both experimental groups
treat.peaks <- intersect(treat1.peaks, treat2.peaks, treatn.peaks)
control.peaks <- intersect(control1.peaks, control2.peaks, controln.peaks)
all.peaks <- union(treat.peaks, control.peaks)

##############################

library("csaw")

# specify paired-end BAMs
treat1.pe.bam <- "treat1.sorted.noDups.filt.noMT.bam"
treat2.pe.bam <- "treat2.sorted.noDups.filt.noMT.bam"
treatn.pe.bam <- "treatn.sorted.noDups.filt.noMT.bam"
control1.pe.bam <- "control1.sorted.noDups.filt.noMT.bam"
control2.pe.bam <- "control2.sorted.noDups.filt.noMT.bam"
controln.pe.bam <- "controln.sorted.noDups.filt.noMT.bam"
pe.bams <- c(treat1.pe.bam, treat2.pe.bam, treatn.pe.bam, control1.pe.bam, control2.pe.bam, controln.pe.bam)

##############################
# read hg38 blacklist
blacklist <- read.table("~/ref_genome/hg38.blacklist.bed", sep="\t")
blacklist <- GRanges(seqnames = blacklist$V1,
                     ranges = IRanges(start = blacklist$V2,
                                      end = blacklist$V3))
# define read parameters
standard.chr <- paste0("chr", c(1:22, "X")) # only use standard chromosomes (and no chrY)
param <- readParam(max.frag=1000, pe="both", discard=blacklist, restrict=standard.chr)

##############################
# get fragment size distribution
treat1.sizes <- getPESizes(treat1.pe.bam)
treat2.sizes <- getPESizes(treat2.pe.bam)
treatn.sizes <- getPESizes(treatn.pe.bam)
control1.sizes <- getPESizes(control1.pe.bam)
control2.sizes <- getPESizes(control2.pe.bam)
controln.sizes <- getPESizes(controln.pe.bam)

# estimate average fragment length by cross-correlation
max.delay <- 1000 # due to 1 kb size-selection step
treat1.corr <- correlateReads(treat1.pe.bam, max.delay, param=param)
treat2.corr <- correlateReads(treat2.pe.bam, max.delay, param=param)
treatn.corr <- correlateReads(treatn.pe.bam, max.delay, param=param)
control1.corr <- correlateReads(control1.pe.bam, max.delay, param=param)
control2.corr <- correlateReads(control2.pe.bam, max.delay, param=param)
controln.corr <- correlateReads(controln.pe.bam, max.delay, param=param)

##############################
# count reads in windows specified by MACS2                                      
peak.counts <- regionCounts(pe.bams, all.peaks, param=param)

##############################
# MACS2 peaks only: filter low abundance peaks
library("edgeR")
peak.abundances <- aveLogCPM(asDGEList(peak.counts))
keep.simple <- peak.abundances > -3 # only use peaks logCPM > -3
peak.counts.filt <- peak.counts[keep.simple,]
# result: few or no peaks should be removed; can modify as desired

##############################
# count BAM reads into 300 bp windows
counts <- windowCounts(pe.bams, width=300, param=param)

# filter uninteresting features (windows) by local enrichment
# local background estimator: neighborhood
surrounds <- 2000
neighbor <- suppressWarnings(resize(rowRanges(counts), surrounds, fix="center"))
wider <- regionCounts(pe.bams, regions=neighbor, param=param)

filter.stat <- filterWindows(counts, wider, type="local")
keep.local <- filter.stat$filter > log2(3) # threshold of 3-fold increase in enrichment over neighborhood abundance; change as desired
counts.local.filt <- counts[keep.local,]

###############################
# count BAM background bins (for TMM normalization)
binned <- windowCounts(pe.bams, bin=TRUE, width=10000, param=param)

##########################################
# NORMALIZATION

# method 1: MACS2 peaks only, TMM normalization based on binned counts
peak.counts.tmm <- peak.counts.filt
peak.counts.tmm <- normOffsets(binned, se.out=peak.counts.tmm)

# method 2: MACS2 peaks only, csaw lowess-normalization
peak.counts.loess <- peak.counts.filt
peak.counts.loess <- normOffsets(peak.counts.loess, type="loess", se.out=TRUE)
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

# method 3: csaw de novo peaks by local enrichment, TMM normalization based on binned counts
counts.local.tmm <- counts.local.filt
counts.local.tmm <- normOffsets(binned, se.out=counts.local.tmm)

# method 4: csaw de novo peaks by local enrichment, csaw lowess-normalization
counts.local.loess <- counts.local.filt
counts.local.loess <- normOffsets(counts.local.loess, type="loess", se.out=TRUE)
# from vignette: "For type="loess", a numeric matrix of the same dimensions as counts, containing the log-based offsets for use in GLM fitting."

#########################################
# DIFFERENTIAL-ACCESSIBILITY ANALYSIS

# set working windows for the desired analysis
working.windows <- peak.counts.tmm # MACS2 peaks only, standard TMM normalization based on binned counts
# working.windows <- peak.counts.loess # MACS2 peaks only, for extreme cases with trended biases (due to global changes in accessibility)
# working.windows <- counts.local.tmm # csaw de novo peaks by local enrichment, standard TMM normalization based on binned counts
# working.windows <- counts.local.loess # csaw de novo peaks by local enrichment, for extreme cases with trended biases (due to global changes in accessibility)
# SEE THE MANUAL FOR MORE INFO ON NORMALIZATION METHODS
###########

# setting up design matrix
# see edgeR manual for more information
y <- asDGEList(working.windows)
colnames(y$counts) <- c("treat1", "treat2", "treatn", "control1", "control2", "controln")
rownames(y$samples) <- c("treat1", "treat2", "treatn", "control1", "control2", "controln")
y$samples$group <- c("treat", "treat", "control", "control")
design <- model.matrix(~0+group, data=y$samples)
colnames(design) <- c("treat", "control")
# design

# stabilise estimates with empirical bayes
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust=TRUE)

# testing for differentially-accessible windows
contrast <- makeContrasts(treat-control, levels=design)
results <- glmQLFTest(fit, contrast=contrast)
# head(results$table)

# combine GRanges rowdata with DE statistics
rowData(working.windows) <- cbind(rowData(working.windows), results$table)
working.windows@rowRanges

# merge nearby windows
# up to "tol" distance apart: 500 bp in this case
# max merged window width: 5000 bp
merged.peaks <- mergeWindows(rowRanges(working.windows), tol=500L, max.width=5000L)
# summary(width(merged.peaks$region))
# result: merges about 300 peaks

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

write.table(final.merged.peaks, "treat_vs_control_macs2_csaw_differential-windows_all.txt",
            sep="\t", quote=F, col.names=T, row.names=F)
write.table(final.merged.peaks.sig, "treat_vs_control_macs2_csaw_differential-windows_significant.txt",
            sep="\t", quote=F, col.names=T, row.names=F)
