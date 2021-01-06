# Sequencing

A set of function to help work with sequencing data (bed files, bam files, fastq files etc...)

## Contains

fromGTF2BED: transforms a GTF file to a BED file, only works for some GTFs for now
getBamDate: parses a bam file header to try to compute when it was generated (as best as it can, if it has had many modification done to it across a long span of time, you will receive the average of that)

## Other very recommended tools

_I am not building anything that overlaps with these tools_

- Bedtools
- samtools
- pyBedtools
- pysam