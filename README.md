# ReFACTor v1.0

Reference-Free Adjustment for Cell-Type composition (ReFACTor) is an unsupervised method for the correction of cell-type heterogeneity in epigenome-wide association studies (EWAS), which is based on a variant of principal component analysis (PCA). ReFACTor is described in the following [paper](http://) (upcoming).

As decribed bellow, ReFACTor is available in R and python. For users working with large datasets we recommend using the much faster python version.

### Download

1. Press the 'Download ZIP' button on the right side of this page
2. Extract the ZIP file to a folder

Dependencies for the python version are desribed at the end of this README file.

### Input

ReFACTor gets the following arguments:

Required:
  * datafile: path to a sites by samples matrix file of beta-normalized methylation levels; the first row should contain the sample IDs and the first column should contain the CpG IDs (see demo_files/demo_datafile.txt for example). Important data preparation  instructions are described below under 'Data preparation'.
  * k: the number of assumed cell types; guidlines for selecting k are desribed below under 'Parameters selection'.

Optional:
  * t: the number of sites to use for computing the ReFACtor components (default is 500); guidlines for selecting k are desribed below under 'Parameters selection'.
  * numcomp: the number of ReFACTor components to output (default is same as k)
  * out: prefix for the output files (default is 'refactor.')

### Output

The software outputs two files:

1. refactor.out.components.txt - a matrix with the first numcomp ReFACtor components for each individual
2. refactor.out.rankedlist.txt - a ranked list of the methylation sites; from the most informative to the least informative

### R

The refactor.R function in the "R" folder implements ReFACTor and can be executed directly from R. For example:

```
# <R code>
source("refactor.R")
k = 5
datafile = "../demo_files/demo_datafile.txt"
results <- refactor(datafile,k)
RC <- results$refactor_components # Extract the ReFACTor components
ranked_list <- results$ranked_list # Extract the list of sites ranked by ReFACTor
# Can also provide optional arguments
results <- refactor(datafile,k,t=500,numcomp=10)
```

##### Demo
The following demo computes the ReFACTor components of a simulated example dataset and performs an EWAS. From the command line run:

```
Rscript demo.R
```


### Python

The refactor.py function in the "python" folder implements ReFACTor and can be executed from the command line as follows:

```
python refactor.py --datafile <datafile> --k <k>
```
or, if including the optional parameters:
```
python refactor.py --datafile <data_file> --k <k> --t<t> --numcomp <numcomp> --out <out_prefix>
```

##### Demo

The following demo computes the ReFACTor components of a simulated example dataset and performs an EWAS. From the command line run:

```
python demo.py
```

### Parameters selection

### Data preparation

##### Preprocessing raw data
ReFACTor expects to get Beta normalized methylation levels (although it may perform well on M-value normalized data as well). Prior to running ReFACTor, the data should be adjusted for known technical atrifacts of the technology used for probing the methylation levels as well as adjusted for known batches. For a comperhensive comparison between methods for preprocessing raw data of the Illumina 450K array see Lenhe et al. 2015, Genome Biology. In order to fit best to the assumptions of ReFACTor, any normalization applied should keep the data approximately normal.

##### Preparing data for ReFACTor
For best performance, we suggest to remove non-autosomal probes, cross-hybridized probes and probes with SNPs. In addition, a large number of sites in the 450K platform are constant or nearly-constant. Removing these sites of very low variance improves the performance of ReFACTor.

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

For any question and for reporting bugs please send an email to Elior Rahmani at: elior.rahmani@gmail.com

