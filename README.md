# ReFACTor v1.0

Reference-Free Adjustment for Cell-Type composition (ReFACTor) is an unsupervised method for the correction of cell-type heterogeneity in epigenome-wide association studies (EWAS), which is based on a variant of principal component analysis (PCA). ReFACTor is described in the following [paper](http://) (upcoming).
s
As decribed bellow, ReFACTor is available in R and in python. For users working with large datasets we recommend using the much faster python version.

### Download

1. Press the 'Download ZIP' button on the right side of the screen
2. Extract the ZIP file to folder

Dependencies for the python version are desribed at the end of this file.

### Input

ReFACTor takes the following arguments:
  * datafile: A data file of sites by samples matrix of beta-normalized methylation levels. The first row should contain the sample IDs and the first column should contain the CpG IDs (see demo_files/demo_datafile.txt for example)
  * k: the number of assumed cell types
  * t (optional): The number of sites to use for computing the ReFACtor components (default is 500)
  * numcomp (optional): The number of ReFACTor components to output (default is same as k)
  * out (optional): Prefix for the output files (default is 'refactor.')

### Output

The software outputs two files:

1. refactor.out.components.txt - a matrix with the first several ReFACtor components for each individual
2. refactor.out.rankedlist.txt - a ranked list of the methylation sites; from the most informative to the least informative

### Execution

#### R

Alternatively, the refactor.R function can be executed within an R script as follows:

```
# <R code>
source("refactor.R")
K = 5
datafile = "demo/demo_datafile.txt"
results <- refactor(datafile,K)
RC <- results$refactor_components # Extract the ReFACTor components
ranked_list <- results$ranked_list # Extract the list of sites ranked by ReFACTor
```

#### python

Execute from the command line:
```
python refactor.py --datafile <datafile> -k <k>
```
or, if including the optional parameters:
```
python refactor.py --datafile <data_file> -k <k> -t<t> -numcomp <num_components> --out <out_prefix>
```

### Demo

The following demo computes the ReFACTor components of a simulated example dataset and performs an EWAS.

#### R

Execute directly from the command line:
```
Rscript demo.R
```

#### python

Execute from the command line:
```
python demo.py
```

### Parameters selection

### Data preprocessing

### Dependencies

For the python version we recommend using a standard python distribution such as Anaconda (https://www.continuum.io/downloads). This release of ReFACTor was implemented for python v2.7 and has the following dependencies:

    numpy
    scipy
    sklearn
    matplotlib (required only for the demo.py)
    statsmodels (required only for the demo.py)

### Citing ReFACTor

If you use ReFACTor in any published work, please cite the manuscript describing the method (upcoming).

### Authors

This software was developed by Reut Yedidim, Noah Zaitlen and Elior Rahmani.

For reporting bugs and questions please email to Elior Rahmani at: elior.rahmani@gmail.com

