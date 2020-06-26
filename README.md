# ATAC-seq

These scripts correspond to a (differential) ATAC-seq analysis workflow as described in [our recent report](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00342-y). It is largely based on the pipeline developed by [Anshul Kundaje's group (Stanford) and the ENCODE project](https://www.encodeproject.org/pipelines/ENCPL792NWO/).

If you use this methodology, please cite the following paper along with corresponding pipeline dependencies:

Jake J. Reske, Mike R. Wilson, and Ronald L. Chandler. 2020. ATAC-seq normalization method can significantly affect differential accessibility analysis and interpretation. *Epigenetics & Chromatin* **13**: 22.

Attempt to run each command individually or in blocks after editing to match your own data architecture.

### ATACseq_workflow.txt
**Generalized ATAC-seq data processing workflow intended for comparative analysis.** Stepwise bioinformatics process and example commands for analyzing ATAC-seq data from raw reads to calling peaks for downstream differential accessibility analysis. Consider “treat1” as an example mouse ATAC-seq Illumina paired-end library. Blue text denotes optional or conditional steps dependent on experimental design and desired output. Users seeking only to discover replicate-concordant accessible regions in a singular cell state may wish to call naïve overlapping peaks, though this step is not necessary for differential accessibility analysis.
![foo](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13072-020-00342-y/MediaObjects/13072_2020_342_Fig4_HTML.png)

### csaw_workflow.R
**csaw workflow for multiple differential accessibility analyses in R.** Graphical representation of proposed csaw workflow in R for calculating differential accessibility. Consider an experimental design with n = 2 biological replicates from two conditions: “treat” and “control”
![foo2](https://media.springernature.com/full/springer-static/image/art%3A10.1186%2Fs13072-020-00342-y/MediaObjects/13072_2020_342_Fig6_HTML.png)
