"""

demo.py demonstrates the utility of the ReFACTor algorithm by performing an EWAS on simulated data in which the phenotype is correlated with the underlying cell type composition (but no true association with the phenotype exists).
The script compares between four different approaches:
- A standatd analysis without adjusting for cell type composition
- Adjusted analysis using the true cell proportions
- Adjusted analysis using ReFACTor
- Adjusted analysis using a standard PCA

The results show that an unadjusted analysis results in an inflation, and while a standard PCA cannot entirely eliminate the inflation, ReFACTor provides the same results as given by the analysis using the real cell proportions.

"""

from numpy import zeros, loadtxt, random, log10, column_stack, ones
import matplotlib.pyplot as plot
from scipy import stats
from methylation_data import MethylationData
import refactor_lib
import statsmodels.api as sm


DATA_FILE = '../demo_files/demo_datafile.txt'             # methylation levels file path
PHENO_FILE = '../demo_files/demo_phenotype.txt'           # phenotype file path
K = 5                                            # the number of assumed cell types
T = 500                                          # number of sites to be selected by ReFACTor
NUM_COMPONENTS = 10                              # number of ReFACTor components to output
OUTPUT_PREFIX = "demo_"                           # prefix for output files names
CELL_COMP_FILE = '../demo_files/demo_cellproportions.txt' # cell composition file path


# run experiments
def run():
    # read methylation data
    meth_data = MethylationData(datafile = DATA_FILE)
    
    # Run ReFACTor; this will output
    refactor  = refactor_lib.Refactor(methylation_data = meth_data, 
                                      K = K, 
                                      t = T, 
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
    draw_qqplot(y=y, title='Unadjusted analysis', xtitle='-log10(expected)', ytitle='-log10(observed)')


    # Run an EWAS corrected for the true cell proportions
    print("Adjusted analysis using cell proportions...")
    R = loadtxt(CELL_COMP_FILE, dtype = float)
    y = associations_test(meth_data, pheno, R)
    plot.subplot(222)
    draw_qqplot(y=y, title='Adjusted analysis using cell proportions', xtitle='-log10(expected)', ytitle='-log10(observed)')


    # Run an EWAS corrected for the first K ReFACTor components
    print("Adjusted analysis using ReFACTor...")
    y = associations_test(meth_data, pheno, refactor.components[:,:K])
    plot.subplot(223)
    draw_qqplot(y=y, title='Adjusted analysis using ReFACTor', xtitle='-log10(expected)', ytitle='-log10(observed)')


    # Run an EWAS corrected for the first PCs of a standard PCA
    print("Adjusted analysis using PCA...")
    y = associations_test(meth_data, pheno, refactor.standard_pca)
    plot.subplot(224)
    draw_qqplot(y=y, title='Adjusted analysis using PCA', xtitle='-log10(expected)', ytitle='-log10(observed)')

    plot.show()

# Generates a QQ-plot for a given vector of p-values.
def draw_qqplot(y, title, xtitle, ytitle, style = 'b.'):
    # x
    unif = random.uniform(size=len(y))
    x = -log10(unif)
    x.sort()

    # y
    y = -log10(y)
    y.sort()
    plot.plot(x, y, style)

    # add y=x trend
    xlim = int(y.max()) + 2
    plot.plot(range(xlim), range(xlim), 'r-') 

    # axis limit
    plot.xlim(0, xlim)
    plot.ylim(0, round(y.max()))

    # titles
    plot.xlabel(xtitle)
    plot.ylabel(ytitle)
    plot.title(title)

# Performs an association test for each methylation site under a linear model.
def associations_test(meth_data, y, model_append = None):
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
    
    return observed_pvalues


if __name__ == "__main__":
    run()

