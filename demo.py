"""

demo.py demonstrates the utility of the ReFACTor algorithm by performing an EWAS on simulated data in which the phenotype is correlated with the underlying cell type composition (but no true association with the phenotype exists).
The script compares between four different approaches:
- A standatd analysis without adjusting for cell type composition
- Adjusted analysis using the true cell proportions
- Adjusted analysis using ReFACTor
- Adjusted analysis using a standard PCA

The results show that an unadjusted analysis results in an inflation, and while a standard PCA cannot entirely eliminate the inflation, ReFACTor provides the same results as given by the analysis using the real cell proportions.

"""

# TODO make sure we explain what this script is doing
# TODO make sure we use all of these imports; same goes for the other files.

from numpy import zeros, loadtxt, random, log10, column_stack, ones
import matplotlib.pyplot as plot
from scipy import stats
from methylation_data import MethylationData
import refactor_lib
import argparse
import statsmodels.api as sm

# Methylation levels
DATA_FILE = 'demo/demo_datafile.txt'
# Phenotype
PHENO_FILE = 'demo/demo_phenotype.txt'
# The number of assumed cell types
K = 5 # TODO const or argument? Elior: keep it const
# Number of sites to be selected by ReFACTor
T = 500
# Number of ReFACTor components to output
NUM_COMPONENTS = 10

OUTPUT_PREFIX = "demo"
CELL_COMP = 'demo/demo_cellproportions.txt'

# run experiment
def run():
    # read methylation data
    meth_data = MethylationData(datafile = DATA_FILE)
    
    # Run ReFACTor; this will output
    refactor  = refactor_lib.Refactor(methylation_data = meth_data, 
                                      K = K, 
                                      t = t, 
                                      num_components = NUM_COMPONENTS, 
                                      ranked_output_filename = OUTPUT_PREFIX + refactor_lib.RANKED_FILENAME, 
                                      components_output_filename  = OUTPUT_PREFIX + refactor_lib.COMPONENTS_FILENAME)

    # Read the phenotype file
    pheno = loadtxt(PHENO_FILE, dtype = str)[:,1:].astype(float)
    pheno = pheno.reshape((meth_data.samples_size,))
   
    # Run an uncorrected EWAS
    print("Unadjusted analysis...")
    y = associations_test(meth_data, pheno)
    plot.subplot(221)
    draw_plot(y=y, style='b.',title='Unadjusted analysis', xtitle='-log10(expected)', ytitle='-log10(observed)', xlim=(0,round(x.max())), ylim=(0,round(y.max())))


    # Run an EWAS corrected for the true cell proportions
    print("Adjusted analysis using cell proportions...")
    R = loadtxt(CELL_COMP, dtype = float)
    y = associations_test(meth_data, pheno, R)
    plot.subplot(222)
    draw_plot(y=y, style='b.',title='Adjusted analysis using cell proportions', xtitle='-log10(expected)', ytitle='-log10(observed)', xlim=(0,round(x.max())), ylim=(0,round(y.max())))


    # Run an EWAS corrected for the first K ReFACTor components
    print("Adjusted analysis using ReFACTor...")
    y = associations_test(meth_data, pheno, self.refactor.components[:,:K])
    plot.subplot(223)
    draw_plot(y=y, style='b.',title='Adjusted analysis using ReFACTor', xtitle='-log10(expected)', ytitle='-log10(observed)', xlim=(0,round(x.max())), ylim=(0,round(y.max())))


    # Run an EWAS corrected for the first PCs of a standard PCA
    print("Adjusted analysis using PCA...")
    #TODO rename the first_pca field to standard_pca
    y = associations_test(meth_data, pheno, refactor.first_pca)
    plot.subplot(224)
    draw_plot(y=y, style='b.', title='Adjusted analysis using PCA', xtitle='-log10(expected)', ytitle='-log10(observed)', xlim=(0,round(x.max())), ylim=(0,round(y.max())))

    plot.show()

# Generates a QQ-plot for a given vector of p-values.
def draw_plot(y, style, title, xtitle, ytitle, xlim, ylim, x=None, line=True):
    if x is None:
        unif = random.uniform(size=self.meth_data.sites_size)
        x = -log10(unif)
        x.sort()
    plot.plot(x, y, 'b.')

    if line: # add y=x line
        plot.plot(range(10), range(10), 'r-') 

    plot.xlabel(xtitle)
    plot.ylabel(ytitle)
    plot.title(title)
    plot.xlim(xlim[0], xlim[1])
    plot.ylim(ylim[0], ylim[1])

# Performs an association test for each methylation site under a linear model.
def associations_test(met_data, y, model_append = None):
    if model_append is not None:
        model_append = column_stack((ones(len(y)), model_append))
    else:
        model_append = ones(len(y))
    observed_pvalues = zeros(meth_data.sites_size)
    for site_i in xrange(meth_data.sites_size):
        c = column_stack((meth_data.data[site_i,:], model_append))
        mod = sm.OLS(y, c)
        res = mod.fit()
        observed_pvalues[site_i] = res.pvalues[0]
    
    #TODO return the p-values and not the -log of the sorted pvalues - the -log and the sort should be done in the draw_plot, which, btw, should be named 'draw_qqplot'.
    y = -log10(observed_pvalues)
    y.sort()
    return y


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', required = True, type = int, help = "The number of assumed cell types")

    args = parser.parse_args()
    run()

