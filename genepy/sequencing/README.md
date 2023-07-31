# Sequencing

A set of function to help work with sequencing data (bed files, bam files, fastq files etc...)

## Contains

- fromGTF2BED: transforms a GTF file to a BED file, only works for some GTFs for now
getBamDate: parses a bam file header to try to compute when it was generated (as best as it can, if it has had many modification done to it across a long span of time, you will receive the average of that)
- getBamDate: from bam files (could be in a google bucket) returns their likely sequencing date if available in the header
- indexBams: given a bucket path, will index all .bam files without an associated index and return their paths
- dropWeirdChromosomes: given a bedfile path, removes chromosomes that are not one of chroms
- extractPairedSingleEndFrom: given a folder, find fastq files and sorts paired and single end based on the R1/R2 patterns
- findReplicates: creates a dict of name and replicate files given a regexp naming scheme
- mergeBams: uses samtools to merge a set of replicates considered into one file

## Other very recommended tools

_I am not building anything that overlaps with these tools_

- Bedtools
- samtools
- pyBedtools
- pysam