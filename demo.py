from numpy import zeros, loadtxt, random, log10, column_stack, ones
import matplotlib.pyplot as plot
from scipy import stats
from methylation_data import MethylationData
import refactor_lib
import argparse
import statsmodels.api as sm

# arguments
DATA_FILE = 'demo/demo_datafile.txt',
K = 5 # TODO const or argument?
T = 500,
PHENO_FILE = 'demo/demo_phenotype.txt',
CELL_COMP = 'demo/demo_cellproportions.txt'):
NUM_COMPONENTS = None, 
OUTPUT_PREFIX = "",

# run experiment
def run():
    # read methylation data
    meth_data = MethylationData(datafile = DATA_FILE)
    
    # run refactor
    refactor  = refactor_lib.Refactor(methylation_data = meth_data, 
                                      K = K, 
                                      t = t, 
                                      num_components = NUM_COMPONENTS, 
                                      ranked_output_filename = output_prefix + refactor_lib.RANKED_FILENAME, 
                                      components_output_filename  = output_prefix + refactor_lib.COMPONENTS_FILENAME)

    # read phenotype file
    pheno = loadtxt(PHENO_FILE, dtype = str)[:,1:].astype(float)
    pheno = pheno.reshape((meth_data.samples_size,))
   
    # exp1
    print("START EXP1")
    y = associations_test(meth_data, pheno)
    plot.subplot(221)
    draw_plot(y=y, style='b.',title='original', xtitle='uniform distribution', ytitle='observed pvalues', xlim=(0,round(x.max())), ylim=(0,round(y.max())))


    # exp2
    print("START EXP2")
    R = loadtxt(CELL_COMP, dtype = float)
    y = associations_test(meth_data, pheno, R)
    plot.subplot(222)
    draw_plot(y=y, style='b.',title='R', xtitle='uniform distribution', ytitle='observed pvalues', xlim=(0,round(x.max())), ylim=(0,round(y.max())))


    # exp3
    print("START EXP3")
    y = associations_test(meth_data, pheno, self.refactor.components[:,:K])
    plot.subplot(223)
    draw_plot(y=y, style='b.',title='refactor', xtitle='uniform distribution', ytitle='observed pvalues', xlim=(0,round(x.max())), ylim=(0,round(y.max())))


    # exp4
    print("START EXP4")
    y = associations_test(meth_data, pheno, refactor.first_pca)
    plot.subplot(224)
    draw_plot(y=y, style='b.', title='P', xtitle='uniform distribution', ytitle='observed pvalues', xlim=(0,round(x.max())), ylim=(0,round(y.max())))

    plot.show()

def draw_plot(y, style, title, xtitle, ytitle, xlim, ylim, x=None, line=True):
    if x is None:
        unif = random.uniform(size=self.meth_data.sites_size)
        x = -log10(unif)
        x.sort()
    plot.plot(x, y, 'b.')

    if line: # add x=y line
        plot.plot(range(10), range(10), 'r-') 

    plot.xlabel(xtitle)
    plot.ylabel(ytitle)
    plot.title(title)
    plot.xlim(xlim[0], xlim[1])
    plot.ylim(ylim[0], ylim[1])

def associations_test(met_data, y, model_append = None):
    if model_append is not None:
        model_append = column_stack((ones(len(y)), model_append))
    else:
        model_append = ones(len(y))

    # get p_value for each site
    observed_pvalues = zeros(meth_data.sites_size)
    for site_i in xrange(meth_data.sites_size):
        c = column_stack((meth_data.data[site_i,:], model_append))

        # linear regression
        mod = sm.OLS(y, c)
        res = mod.fit()
        
        # find p_values    
        observed_pvalues[site_i] = res.pvalues[0]
    
    
    y = -log10(observed_pvalues)
    y.sort()
    return y


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', required = True, type = int, help = "The number of assumed cell types")

    args = parser.parse_args()
    run()

