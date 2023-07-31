# cell_line_mapping
Code for mapping between different CCLE/DepMap cell line identifiers

## Installation
### R
```
options(repos = c(
	"https://iwww.broadinstitute.org/~datasci/R-packages",
	"https://cran.cnr.berkeley.edu"))
install.packages('celllinemapr')
```
As one could not have access to intranet to download the name mapping. The name mapping file is directly available and can be put to work by executing this command:
`mkdir ~/.celllinemapr && mkdir ~/.celllinemapr/data && cp naming.csv ~/.celllinemapr/data`

### Python
```
pip install https://intranet.broadinstitute.org/~datasci/python-packages/cell_line_mapper-latest.tar.gz
```

## Usage
See [here](https://github.com/broadinstitute/cell_line_mapping/blob/master/celllinemapr/SOP.Rmd) for examples of functions for the R package.
Funtion names for the Python package are analogous to those for R, replacing `.` with `_`. For instance,

R: `ccle.to.arxspan`

Python: `ccle_to_arxspan`
