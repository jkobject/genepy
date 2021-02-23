# JKBio

A set of awesome functions & tools for Computational Genomists

![long genome](documentation/genome.jpg)
## Content

- **utils**: where a bunch of helper functions and usefull general scripts are stoed
  - **plots**: a set of plotting tools based on [matplotlib]() and [bokeh]() to make volcano plots / CNV maps etc..
  - **helper**: and additional helper functions to save data, do merging of dataframes...
- **terra**: contains a set of functions that uses [dalmatian]() to interact with the [GCP]() powered genomics HPC platform: [Terra](). 
- **sequencing**: contains a set of function to works with bed/bam/fastqs...
- **rna**: contains function to work with RNAseq (and related) data.
  - **pyDESeq2**: it is a python integration of [deseq2]() (the differential expression analyser) with [rpy2]()
- **mutations**: a set of functions to work with maf files, vcf files etc..
- **google**: functions and packages linked to google's apis
  - **google_sheet**: function to upload a df as a google sheet
  - **gcp**: sets of functions to interact with google storage (relies on gsutil)
- **epigenetics**: where we have things related to epigenomics
  - **rose**: where an updated version of the rose algorithm is stored (as a git submodule)
  - **chipseq**: has functions to read, merge, denoise, ChIP seq data, it contains a lot of functions required for the AML paper.
- **taigr**: a version of taiga that do not requires RCurl (and can save you when you have a faulty RCurl-Curl link)
- **data**: should not contain anything when pulled but is used by any of the functions in the folder, to save some required data files
- **cell_line_mapping**: a set of functions to map cell line ids to other cell line ids based on an up to date google spreadsheet. 


## Install

### with pip (WIP)

`pip install JKBio`
### dev mode (better for now)

```bash
git clone git://github.com/jkobject/JKBio.git
cd JKBio
git submodule update --init
```

then you can import files in python with e.g:
```python
from JKBio import TerraFunction as terra
```

if JKBio is not in your path, first do:

```python
import sys
sys.path.append(RELATIVE_PATH_TO_JKBio)
```

now you can install the necessary python packages:

```bash
pip install requirements.txt
pip install rpy2-bioconductor-extensions gseapy macs2 deeptools
```

or if not using the requirements.txt (computation results might change):

```bash
pip install numpy pandas
```

```bash
pip install bokeh dalmatian firecloud_dalmatian google_api_python_client gsheets gspread ipdb ipython matplotlib Pillow pybedtools pyBigWig pysam pytest requests rpy2 scikit_learn scipy seaborn setuptools taigapy taigapy typing venn rpy2-bioconductor-extensions gseapy macs2 deeptools
```

then install the following tools:
- [htslib/samtools](http://www.htslib.org/)
- [bwa](https://github.com/lh3/bwa)
just used once:
- [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)

finaly you can install R packages (GSEABase, erccdashboard, GSVA, DESeq2):

```bash
R -e 'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")};BiocManager::install(c("GSEABase", "erccdashboard", "GSVA", "DESeq2"));'
```
## About

As I am working in different domains of computational genomics, I need to have a set of reusable function that will help me during my work.
It can be functions to work with different tools that I have to use. Functions to do some plots. etc..

I will be trying to keep each of these functions functional and documented. Feel free to pull and start use anything that might be useful to you.
If you see anything suspicious or not working. A pull request would definitely get reviewed within a day.

I hope to be able to give back to the community and maybe save a couple of hours to couple of researchers.

Best.


jkalfon@broadinstitute.org

jkobject@gmail.com

https://jkobject.com

Apache license 2.0.
