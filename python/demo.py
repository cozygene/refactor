"""

demo.py demonstrates the utility of the ReFACTor algorithm by performing an EWAS on simulated data in which the phenotype is correlated with the underlying cell type composition (but no true association with the phenotype exists).

The script first runs ReFACTor on the data and outputs two files:
- demo_refactor.out.components.txt: the ReFACTor components
- demo_refactor.out.rankedlist.txt:  a ranked list of the CpGs, from the most informative site to the least informative site (as determined by ReFACTor)

Second, the script performes an EWAS. Four different approaches are compared:
- A standatd analysis without adjusting for cell type composition
- Adjusted analysis using the true cell proportions of the simulated data
- Adjusted analysis using ReFACTor
- Adjusted analysis using a standard PCA

The results of the EWAS are saved into a file (demo_results.png), showing that unadjusted analysis results in an inflation, and while a standard PCA cannot entirely eliminate this inflation, ReFACTor can and provides the same correction given by the analysis that used the real cell proportions.

"""

from numpy import zeros, loadtxt, random, log10, column_stack, ones, linspace
import matplotlib.pyplot as plot
from scipy import stats
from refactor_lib import methylation_data
from refactor_lib import refactor
import statsmodels.api as sm


K = 5                                                     # the number of assumed cell types
# Simulated data:
DATA_FILE = '../demo_files/demo_datafile.txt'             # methylation levels file path
PHENO_FILE = '../demo_files/demo_phenotype.txt'           # phenotype file path
CELL_COMP_FILE = '../demo_files/demo_cellproportions.txt' # cell composition file path

def run():

    # Run ReFACTor
    refactor_obj  = refactor.Refactor(DATA_FILE, K, out="demo_refactor")

    # read methylation data
    meth_data = methylation_data.MethylationData(DATA_FILE)

    # Read the phenotype file
    pheno = loadtxt(PHENO_FILE, dtype = str)[:,1:].astype(float)
    pheno = pheno.reshape((meth_data.samples_size,))
    
    # set plots
    fig, axes = plot.subplots(nrows=4, ncols=4)
    fig.set_tight_layout(True)
    
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


    # Run an EWAS corrected for the first k ReFACTor components
    print("Adjusted analysis using ReFACTor...")
    y = associations_test(meth_data, pheno, refactor_obj.components[:,:K])
    plot.subplot(223)
    draw_qqplot(y=y, title='Adjusted analysis using ReFACTor', xtitle='-log10(expected)', ytitle='-log10(observed)')


    # Run an EWAS corrected for the first PCs of a standard PCA
    print("Adjusted analysis using PCA...")
    y = associations_test(meth_data, pheno, refactor_obj.standard_pca)
    plot.subplot(224)
    draw_qqplot(y=y, title='Adjusted analysis using PCA', xtitle='-log10(expected)', ytitle='-log10(observed)')
    
    plot.savefig("demo_results.png")
    print("Plotted and saved the results into demo_results.png")

    # plot.show()

# Generates a QQ-plot for a given vector of p-values.
def draw_qqplot(y, title, xtitle, ytitle, style = 'b.'):
    # x
    #unif = random.uniform(size=len(y))
    #x = -log10(unif)
    #x.sort()
    x = -log10(linspace(0,1,len(y)+1))
    x = x[1:]
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

