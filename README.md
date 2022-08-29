# genepy

_what is [genepy](https://en.wikipedia.org/wiki/G%C3%A9n%C3%A9pi)?_

A set of awesome functions & tools for Computational Geneticists

![long genome](documentation/genome.jpg)

## Content

- **utils**: where a bunch of helper functions and usefull general scripts are stored
  - **plots**: a set of plotting tools based on [matplotlib]() and [bokeh]() to make volcano plots / CNV maps etc..
  - **helper**: and additional helper functions to save data, do merging of dataframes...
- **terra**: contains a set of functions that uses [dalmatian]() to interact with the [GCP]() powered genomics HPC platform: [Terra](). 
- **sequencing**: contains a set of function to works with bed/bam/fastqs...
- **rna**: contains function to work with RNAseq (and related) data.
  - **pyDESeq2**: it is a python integration of [deseq2]() (the differential expression analyser) with [rpy2]()
- **mutations**: a set of functions to work with maf files, vcf files etc..
- **google**: functions and packages linked to google's apis
  - **google_sheet**: function to upload a df as a google sheet
  - **gcp**: sets of functions to interact with google storage (relies on `gsutil`)
- **epigenetics**: where we have things related to epigenomics
  - **chipseq**: has functions to read, merge, denoise, ChIP seq data.
  - **plot**: has functions to plot ChIP seq data.

### Helper tools

_tools that you do not need to use directly as they have binding functions in genepy._ 

- **epigenetics/rose:**: where an updated version of the rose algorithm is stored (as a git submodule) 
- **cell_line_mapping-master/python/cell_line_mapper**: a set of functions to map cell line ids to other cell line ids based on an up to date google spreadsheet. 


## Install

### with pip

`pip install broad-genepy`

and then use with `from genepy.utils/epigenetics/... import ...`

Please see the next step to get access to all bindings and tools.

### dev mode

```bash
git clone git://github.com/BroadInstitute/genepy.git
pip install -e genepy
```

then you can import files in python with e.g:
```python
from genepy import terra
from genepy.utils import helper as h
from genepy.google import gcp
from genepy.utils import plot
from genepy.epigenetics import chipseq

```

## installation: to get access to all bindings and tools

Install the following tools:
- [gcloud](https://cloud.google.com/sdk/docs/install-sdk)
- [firecloud-dalmatian](https://github.com/getzlab/dalmatian) 
- [gsheets](https://github.com/xflr6/gsheets)
- [htslib/samtools](http://www.htslib.org/)
- [bwa](https://github.com/lh3/bwa)
just used once:
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

Some of these packages like gsheets, gcloud, firecloud-dalmatian will require you to create google accounts, login on your machine or download oauth files.

Finaly you can install R packages (GSEABase, erccdashboard, GSVA, DESeq2):

```bash
R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};BiocManager::install(c("GSEABase", "erccdashboard", "GSVA", "DESeq2"));'
```

## data:

hg38 genome sizes: from https://github.com/igvteam/igv/blob/master/genomes/sizes/hg38.chrom.sizes

## About

please do contribute, we do not have time to fix all issues or work on feature requests

Jeremie Kalfon jkalfon@broadinstitute.org jkobject@gmail.com https://jkobject.com

Apache license 2.0.
