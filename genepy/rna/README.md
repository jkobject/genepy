# RNA

A set of functions to work with RNAseq (and related) data type

## contains


- filterProteinCoding: removes all non protein coding genes from a list (you need taiga access)
- convertGenes: converts genes from a naming to another (you need taiga access)
- getSpikeInControlScales: extracts the spike in control values from a set of bam files
- GSEAonExperiments: perform GSEA to compare a bunch of conditions at once
- runERCC: creates an ERCC dashboard and extract the RNA spike ins from it (need rpy2 and ipython and R's ERCCdashboard installed)

## recommended tools

- ERCCdashboard (R)
- DESeq2 (R)
- slamdunk
- GSVA (R)
- gseapy (python)