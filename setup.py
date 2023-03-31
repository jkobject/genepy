from setuptools import setup
import sys
import os
import io
import subprocess

if sys.version_info.major < 3 or sys.version_info.minor < 2:
    raise ValueError("genepy is only compatible with Python 3.5 and above")
if sys.version_info.minor < 5:
    import warnings

    warnings.warn("genepy may not function properly on Python < 3.8")

os.system("git submodule init && git submodule sync")

with open("README.md", "r") as f:
    long_description = f.read()

print("trying to install R packages")
try:
    subprocess.run(
        'R -e \'if(!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager", repos="http://cran.us.r-project.org")};BiocManager::install(c("GSEABase", "erccdashboard", "GSVA", "DESeq2"));\'',
        shell=True,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    subprocess.run("pip install rpy2")
except:
    print("R packages not installed")
print("if it did not work. please install R or check your R installation")
print(
    "once R is installed you need to install erccdashboard, GSEABase GSVA, DESeq2 to have access to all the functions"
)


def read(*paths, **kwargs):
    """Read the contents of a text file safely.
    >>> read("genepy", "VERSION")
    '0.1.0'
    >>> read("README.md")
    ...
    """

    content = ""
    with io.open(
        os.path.join(os.path.dirname(__file__), *paths),
        encoding=kwargs.get("encoding", "utf8"),
    ) as open_file:
        content = open_file.read().strip()
    return content


def read_requirements(path):
    return [
        line.strip()
        for line in read(path).split("\n")
        if not line.startswith(('"', "#", "-", "git+"))
    ]


setup(
    name="wlg-genepy",
    version=read("genepy", "VERSION"),
    description="A useful module for any CompBio",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Jeremie Kalfon",
    author_email="jkobject@gmail.com",
    url="https://github.com/BroadInstitute/genepy",
    packages=[
        "genepy/cell_line_mapping-master/python/cell_line_mapper",
        "genepy/epigenetics",
        "genepy/mutations",
        "genepy/google",
        "genepy/sequencing/",
        "genepy/terra",
        "genepy/rna",
        "genepy/utils",
    ],
    package_data={"genepy": ["data/*"]},
    python_requires=">=3.5",
    install_requires=read_requirements("requirements.txt"),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)


print(
    "You might want to install Bowtie2, samtools, bwa and R to be able to use all functions of this package:\n\
  http://bowtie-bio.sourceforge.net/bowtie2/index.shtml\n\
  http://www.htslib.org/\n\
  https://github.com/lh3/bwa\n"
)

print("Finished!")
