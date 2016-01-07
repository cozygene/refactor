from numpy import zeros, loadtxt, random, log10, column_stack
import matplotlib.pyplot as plot
from scipy import stats
# from regression import LinearRegression
from methylation_data import MethylationData
import refactor_lib
import argparse
import os
import sys
import configurelogging
configurelogging.configureLogging('')
import logging
import statsmodels.api as sm

class RefactorDemo(object):

    def __init__(self,
                  K,
                  datafile = 'demo/demo_datafile.txt',
                  t = 500,
                  num_components = None, 
                  output_prefix = "",
                  phenofile = 'demo/demo_phenotype.txt',
                  cellcomp = 'demo/demo_cellproportions.txt'):

        meth_data = MethylationData(datafile = datafile)
        self.refactor  = refactor_lib.Refactor(methylation_data = meth_data, 
                                                K = K, 
                                                t = t, 
                                                num_components = num_components, 
                                                ranked_output_filename = output_prefix + refactor_lib.RANKED_FILENAME, 
                                                components_output_filename  = output_prefix + refactor_lib.COMPONENTS_FILENAME)

        self.meth_data = meth_data
        self.y = self._load_and_validate_phenotype(phenofile).reshape((self.meth_data.samples_size,))
        self.R = loadtxt(cellcomp, dtype = float)

        # exp1
        xasix = range(10)
        yasix = range(10)
        # _, plots_arr = plot.subplots(2, 2)
        print("START EXP1")
        x,y = self.associations_test(self.y, zeros((self.meth_data.samples_size ,)))
        plot.subplot(221)
        plot.plot(x, y, 'b.')
        plot.plot(stats.linregress(xasix, yasix))
        plot.xlabel('uniform distribution')
        plot.ylabel('observed pvalues')
        plot.title('original')
        plot.xlim(0,round(x.max()))
        plot.ylim(0,round(y.max()))
        # plot(res$x, res$y, main="original", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
        # abline(lm(xasix~yasix),col=2,lty=1)


        # exp2
        print("START EXP2")
        x,y = self.associations_test(self.y, self.R)
        plot.subplot(222)
        plot.plot(x, y, 'b.')
        plot.plot(stats.linregress(xasix, yasix))
        plot.xlabel('uniform distribution')
        plot.ylabel('observed pvalues')
        plot.title('R')
        plot.xlim(0,round(x.max()))
        plot.ylim(0,round(y.max()))
        # plot(res$x, res$y, main="R", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
        # abline(lm(xasix~yasix),col=2,lty=1)

        # exp3
        print("START EXP3")
        
        x,y  = self.associations_test(self.y, self.refactor.components[:,:K])
        plot.subplot(223)
        plot.plot(x, y, 'b.')
        plot.plot(stats.linregress(xasix, yasix))
        plot.xlabel('uniform distribution')
        plot.ylabel('observed pvalues')
        plot.title('refactor')
        plot.xlim(0,round(x.max()))
        plot.ylim(0,round(y.max()))
        # plot(res$x, res$y, main="refactor", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
        # abline(lm(xasix~yasix),col=2,lty=1)

        # exp4
        print("START EXP4")
        x,y = self.associations_test(self.y, self.refactor.first_pca)
        plot.subplot(224)
        plot.plot(x, y, 'b.')
        plot.plot(stats.linregress(xasix, yasix))
        plot.xlabel('uniform distribution')
        plot.ylabel('observed pvalues')
        plot.title('P')
        plot.xlim(0,round(x.max()))
        plot.ylim(0,round(y.max()))
        # plot(res$x, res$y, main="P", xlab="uniform distribution", ylab="observed pvalues", pch='.', xlim=c(0,ceiling(max(res$x))),  ylim=c(0,ceiling(max(res$y))))
        # abline(lm(xasix~yasix),col=2,lty=1)

        # dev.off()
        # print("plot saved to plot.png file") #TODO move plot.png to argument
        #                                      #TODO change plot.png by time and date
        # system("open plot.png")

        plot.show()


    def associations_test(self, y, model_append):
        observed_pvalues = zeros((self.meth_data.sites_size,))
        import pdb
        # pdb.set_trace()
        for site_i in xrange(self.meth_data.sites_size):
            c =   column_stack((self.meth_data.data[0,:], model_append))
            # model = LinearRegression(y, c) # find most effocoent func
            mod = sm.OLS(y, c)
            res = mod.fit()
            
            # slope, intercept, r_value, p_value, std_err = stats.linregress(y, c)
            # TODO take only firt p_value
            # pvalue <- coef(summary(model))[2,4]
            observed_pvalues[site_i] = res.pvalues[0]
        
        unif = random.uniform(size=self.meth_data.sites_size)
        # x <- runif(m)
        x = -log10(unif)
        y = -log10(observed_pvalues)
        x.sort()
        y.sort()
        return x, y

    def _load_and_validate_phenotype(self, phenofile):
        pheno = None
        if phenofile:
            pheno = self._load_and_validate_matrix_ids_and_reorder(phenofile)
            if len(pheno[0]) != 2:
                logging.error("must provided only one phenotype. should be 2 columns: 1 - sample id, 2 - phenotype") #TODO is this right?
                sys.exit(2)

            pheno = pheno[:,1:].astype(float) # TODO should check if can convert  to float
            

        return pheno


    """
    reads matrix from matrix_file_path
    validates that the matrix has number of rows as the number of sample ids
    checks that the sample ids in matrix (the first column) are the same ids as in sample_ids list
    if they are the same but in different order, reorder matrix rows as in sample_ids
    """
    def _load_and_validate_matrix_ids_and_reorder(self, matrix_file_path):
        if not os.path.exists(matrix_file_path) :
            logging.error("The file '%s' doesn't exist. Exiting" % matrix_file_path)
            sys.exit(2)

        data = loadtxt(matrix_file_path, dtype = str)
        if len(data) != len(self.meth_data.samples_ids):
            logging.error("the file provided %s doesn't include all sample ids" % matrix_file_path)
            sys.exit(2)

        matrix_sample_ids = data[:,0]

        if not (self.meth_data.samples_ids == matrix_sample_ids).all(): #todo check this is not by order
            if len(set(self.meth_data.samples_ids)^set(matrix_sample_ids)) != 0:
                logging.error("sample ids in phenotype file are not the same as in the data file")
                sys.exit(2)
            
            logging.info("sample ids in phenotype file are not in the same order as in the datafile, reordering...")
            sample_to_row = dict()
            for sample in data:
                sample_to_row[sample[0]] = sample

            orderd_data = empty_like(data)  
            for i,sid in enumerate(self.meth_data.samples_ids):
                orderd_data[i,:] = sample_to_row[sid]
            return orderd_data

        return data

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', required = True, type = int, help = "The number of assumed cell types")

    args = parser.parse_args()

    RefactorDemo(args.k)

