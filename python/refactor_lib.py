import os
import sys
from sklearn import preprocessing
from numpy import dot, linalg, sqrt, loadtxt
import pca
from regress import Regress
import matplotlib.pyplot as plot
from numpy import dot, linalg, sqrt, loadtxt, linalg, log
from methylation_data import MethylationData

class Refactor( object ):
    
    RANKED_FILENAME =       '.out.rankedlist.txt'
    COMPONENTS_FILENAME =   '.out.components.txt'
    VERSION = 1.0 

    def __init__( self,
                  datafile,
                  k,
                  covarfile = None,
                  t = 500,
                  num_components = None,
                  stdth = 0.02,
                  out = "refactor"
                ):

        # load methylation data file
        self.meth_data = MethylationData(datafile)
        
        # validate and process all variables
        self.k =                          self._validate_k(k)
        self.t =                          self._validate_t(t)
        self.num_components =             self._validate_num_comp(num_components)
        self.ranked_output_filename =     out + self.RANKED_FILENAME
        self.components_output_filename = out + self.COMPONENTS_FILENAME

        self.run(covarfile, stdth)


    def _validate_file_path(self, filepath):
        if not os.path.exists(filepath) :
            print("ERROR: The file '%s' doesn't exist. Exiting" % filepath)
            self._terminate_refactor()

    def _validate_k(self,k):
        if not (k >= 2 and k <= self.meth_data.samples_size):
            print("ERROR: k must be at least 2 and smaller than the number of samples. k = %s, samples = %s" % (k, self.meth_data.samples_size))
            self._terminate_refactor()

        return k

    def _validate_t(self,t):
        if t > self.meth_data.sites_size or t < self.k : 
            print("ERROR: t cannot be greater than the number of sites or smaller than k . t = %s, sites = %s, k = %s" % (t, self.meth_data.sites_size, self.k))
            self._terminate_refactor()

        return t

    def _validate_num_comp(self,num_comp):
        if num_comp and not (num_comp >= self.k and num_comp <= self.meth_data.samples_size):
            print("ERROR: the number of components must be at least k and smaller than the number of samples. num_comp = %s, samples = %s, k = %s" % (self.t, self.meth_data.samples_size, self.k))
            self._terminate_refactor()

        return num_comp if num_comp else self.k

    def _load_and_validate_covarfile(self, filepath):
      ids = loadtxt(filepath, dtype = str)[:,0].astype(str)
      if len(ids) != len(self.meth_data.samples_ids) or sum(self.meth_data.samples_ids[i] == ids[i] for i in range(len(ids))) < len(ids):
            print("ERROR: The order of the samples in the covariates file must be the same as the order in the data file")
            self._terminate_refactor()
      covs = loadtxt(filepath, dtype = str)[:,1:].astype(float)
      return covs

    def _adjust_data(self, covs):
        data_adj = self.meth_data._copy()
        for i in range(self.meth_data.sites_size):
            data_adj.data[i,:] = Regress.regress(self.meth_data.data[i,:],covs)
        return data_adj

    def run(self, covarfile, stdth):
        print("Excluding sites with low variance (std < %s)..." % str(stdth))
        before = self.meth_data.sites_size
        self.meth_data._filter_sites_by_std(stdth)
        print("%d sites were excluded due to low variance..." % (before - self.meth_data.sites_size) )

        if (covarfile is not None):
            print('Adjust data for covariates...')
            covs = self._load_and_validate_covarfile(covarfile)
            self.meth_data = self._adjust_data(covs)

        print('Starting ReFACTor v%s...' % self.VERSION);
        self.components, self.ranked_sites, self.standard_pca = self._refactor()
        print('ReFACTor is done!')

   
    def _write_file( self, filepath, data):   
        if  os.path.exists(filepath):
            os.remove(filepath) 
        with open(filepath, 'w') as f:
            f.write(data)

    def _refactor( self ):
        print('Running a standard PCA...')
        pca_out1 = pca.PCA(self.meth_data.data.transpose()) 

        print('Computing a low rank approximation of the input data and ranking sites...')
        x = self._low_rank_approximation(pca_out1.P, pca_out1.U, self.k)
        
        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(self.meth_data.data.transpose()).transform(self.meth_data.data.transpose())
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x).transform(x)
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )

        distances = self._euclidean_distance(An, Bn)
        ranked_list = distances.argsort()

        print('Computing the ReFACTor components...')
        sites = ranked_list[0:self.t]

        pca_out2 = pca.PCA(self.meth_data.data[sites,:].transpose())
        score = pca_out2.P

        print('Saving a ranked list of the data features...')
        #data = '\n'.join(['%s\t%s'% (index+1, self.meth_data.cpgnames[index]) for index in ranked_list])
        data = '\n'.join(['%s'% self.meth_data.cpgnames[index] for index in ranked_list])
        self._write_file(self.ranked_output_filename, data)

        print('Saving the ReFACTor components...')
        data = '\n'.join(['\t'.join([str(i) for i in line]) for line in score[:,0:self.num_components   ]])
        self._write_file(self.components_output_filename, data)
        
        return score[:,0:self.num_components], ranked_list, pca_out1.P[:,0:self.k]


    def _low_rank_approximation(self, A, B, i):
        return dot(A[:,0:i], B[:,0:i].transpose())

    def _euclidean_distance(self, A, B):
        return sqrt(((A - B)**2).sum(axis=0))

    def _terminate_refactor(self):
        print("ReFACTor was terminated.")
        exit(2)


    @staticmethod
    def estimate_k(methylation_data, max_k):

        min_k = 2

        # Find the eigenvalues of the covariance matrix
        eigs = sorted(linalg.eigvals(dot(methylation_data.data.transpose(),methylation_data.data)),key=lambda x: -x)

        # For each eigenvalue i compute its score: -log of the ratio between the i-th eigenvalue and the (i-1)-th eigenvalue.
        scores = [0 for i in range(max_k-min_k+1)]
        counter = 0
        for i in range(min_k,max_k):
            scores[counter] = -log(eigs[i-1] / eigs[i-2])
            counter += 1

        # Plot #eigenvalue vs. scores
        fig, axes = plot.subplots(nrows=1, ncols=1)
        plot.plot([i for i in range(min_k,max_k+1)], scores)
        plot.xlabel('# eigenvalue')
        plot.ylabel('score')
        filename = "estimate_k_results.png"
        plot.savefig(filename)
        print("Plotted and saved the results into %s" % filename)


    @staticmethod
    def estimate_t(methylation_data, k, numsites):

        span = 9 # parameter for the moving average (the window size for constructing the average); must be an even number.

        # Compute a low rank approximation of the data
        pca_res = pca.PCA(methylation_data.data.transpose()) 
        x = dot(pca_res.P[:,0:k], pca_res.U[:,0:k].transpose())        
        
        # Compute the distance of each site form its low rank approximation
        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(methylation_data.data.transpose()).transform(methylation_data.data.transpose())
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x).transform(x)
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )
        distances = sorted(sqrt(((An - Bn)**2).sum(axis=0)))

        # Compute a score for each site i of the sorted distances list: the moving average of dist(i) - dist(i-1)
        distances_diff = [distances[0]] + [distances[i] - distances[i-1] for i in range(1,numsites)]

        scores = [0 for i in range(numsites)]
        mid = (span-1) / 2        
        for i in range(0,numsites):
            l = [distances_diff[j] for j in range(min(0,i-mid),min(i+mid,numsites))]
            scores[i] = sum(l) / float(len(l))
        
        # Plot sites vs. scores
        fig, axes = plot.subplots(nrows=1, ncols=1)
        plot.plot([i for i in range(1,numsites+1)], scores)
        plot.xlabel('site')
        plot.ylabel('score')
        filename = "estimate_t_results.png"
        plot.savefig(filename)
        print("Plotted and saved the results into %s" % filename)
        

