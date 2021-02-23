from setuptools import setup

with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='JKBio',
    version='1.0',
    description='A useful module for any CompBio',
    long_description = long_description,
    author='Jeremie Kalfon',
    author_email='jkobject@gmail.com',
    url="https://github.com/jkobject/JKBio",
    packages=['cell_line_mapping-master/python','epigenetics/rose'],  # same as name
    python_requires='>=3.5',
    install_requires=[
        'rpy2-bioconductor-extensions',
        'gseapy',
        'macs2',
        'deeptools',
        ## from requirements.txt
        "bokeh",
        "dalmatian",
        "firecloud_dalmatian",
        "google_api_python_client",
        "gsheets",
        "gspread",
        "ipdb",
        "ipython",
        "matplotlib",
        "numpy",
        "pandas",
        "Pillow",
        "pybedtools",
        "pyBigWig",
        "pysam",
        "pytest",
        "requests",
        "rpy2",
        "scikit_learn",
        "scipy",
        "seaborn",
        "setuptools",
        "taigapy",
        "taigapy",
        "typing",
        "venn",
        ],  # external packages as dependencies
)

print("You might want to install Bowtie2, samtools, bwa and R to be able to use all functions of this package:\n\
  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml\n\
  http://www.htslib.org/\n\
  https://github.com/lh3/bwa\n")

print("once R is installed you need to have installed erccdashboard, GSEABase GSVA, DESeq2 to have access to aall the functions")

print("Finished!")
