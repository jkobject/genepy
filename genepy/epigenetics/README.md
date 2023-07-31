# epigenomics

Especially targeted to functions related to the analysis of epigenomics data. It has functions to read, merge, denoise, ChIP seq data.

## Available functions:

### chipseq.py

- bigWigFrom: run the bigwig command line for a set of bam files in a folder
- ReadRoseSuperEnhancers: reads ROSE2's output and returns its superenhancer bedfile as a pd dataframe. 
- loadPeaks: loads 1 to many peak bedfile into one pandas dataframe.
- simpleMergePeaks: simply merges bedfiles from peak callers. providing a concaneted dataframe of bed-like tables
- putInBed: given a conscensus bed-like dataframe and another one, will merge the second one into the first
- pairwiseOverlap: compute pairwise overlap and correlation on this overlap for a set of peaks mappe to a conscensus 
- enrichment: compute pairwise enrichment and correlation for a set of peaks mappe to a conscensus 
- fullDiffPeak: will use macs3 to call differential peak binding from two bam files and their control
- diffPeak: calls MACS2 bdgdiff given some parameters
- MakeSuperEnhancers: Calls super enhancer from H3K27ac with the ROSE algorithm
- runChromHMM: runs the chromHMM algorithm
- loadMEMEmotifs: loads motif from the output file of MEME after running fimo.
- simpleMergeMotifs: aggregates the motifs if they overlap, into one motif file
- substractPeaksTo: removes all peaks that are not within a bp distance to a set of loci

### CREME.py

The goal of creme is to be a simple, 1 function tool. It works with 1 to many sets of replicates for each pulled protein/mark.

CREME takes as input a pandas dataframe. This dataframe is the concatenation of each replicates' bed files and can be loaded from a set of MACS2 bedfiles using genepy's loadPeaks function.

CREME will output, amongst other thing, a dataframe representing a concatenation f bedfiles of merged replicates.

find out more at __CREME.md__

## highly recommended packages

*This package won't contain anything that overlap with those and might use those packages for what it is doing.*
- Bedtools
- deepTools
- MACS2
- ROSE
- MEME
- ChromHMM
