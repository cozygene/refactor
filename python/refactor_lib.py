import os
import sys
from sklearn import preprocessing
from numpy import dot, linalg, sqrt, loadtxt
import pca

#TODO NOTE remember to copy the matrix before making changes!!!!

RANKED_FILENAME =       'refactor.out.rankedlist.txt'
COMPONENTS_FILENAME =   'refactor.out.components.txt'

class Refactor( object ):
    VERSION = 1.0 #TODO move to other place?

    def __init__( self,
                  methylation_data,
                  K,
                  t = 500,
                  num_components = None, 
                  ranked_output_filename = RANKED_FILENAME,
                  components_output_filename = COMPONENTS_FILENAME
                ):
        """
        methylation_data is a MethylationData object
        TODO add class doc here
        """
        
        self.meth_data = methylation_data
        
        # validate and process all variables
        self.k =                          self._validate_k(K)
        self.t =                          self._validate_t(t)
        self.num_components =             self._validate_num_comp(num_components)
        self.ranked_output_filename =     ranked_output_filename
        self.components_output_filename = components_output_filename

        self.run()



    def _validate_file_path(self, filepath):
        if not os.path.exists(filepath) :
            logging.error("The file '%s' doesn't exist. Exiting" % filepath)
            sys.exit(2)

    def _load_and_validate_file_of_dimentions(self, filepath, dim):
        """
        validates that a file exists and that it is a matrix from dimentions dim
        """
        if filepath is None:
            return None

        self._validate_file_path(filepath)
        print("Loading file %s..." % filepath)
        data = loadtxt(filepath, dtype = str)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if len(data.shape) != dim:
            logging.error("The file '%s' is not a %sd matrix" % (filepath, dim))
            sys.exit(2)

        return data

    """
    2 <= K <= samples size
    """
    def _validate_k(self,k):
        if not (k >= 2 and k <= self.meth_data.samples_size):
            logging.error("k must be at least 2 and smaller than samples size. k = %s, samples = %s" % (k, self.meth_data.samples_size))
            sys.exit(2) 

        return k

    """
    K <= t <= sites size
    must be called after _validate_k
    """
    def _validate_t(self,t):
        if t > self.meth_data.sites_size or t < self.k : 
            logging.error("t cannot be greater than number of sites or smaller than k . t = %s, sites = %s, k = %s" % (t, self.meth_data.sites_size, self.k))
            sys.exit(2) 

        return t

    """
    K <= num_comp  <= samples size
    must be called after _validate_k
    """
    def _validate_num_comp(self,num_comp):
        if num_comp and not (num_comp >= self.k and num_comp <= self.meth_data.samples_size):
            logging.error("number of components must be at least k and smaller than samples size. num_comp = %s, samples = %s, k = %s" % (t, self.meth_data.samples_size, self.k))
            sys.exit(2) 

        return num_comp if num_comp else self.k

    def run( self ):
        print('Starting ReFACTor v%s...' % self.VERSION);
        self.components, self.ranked_sites, self.standard_pca = self._refactor()
        print('ReFACTor Done!')

   
    """
    writes data to file filepath.
    removes the file if already exists
    """
    def _write_file( self, filepath, data):   
        if  os.path.exists(filepath):
            os.remove(filepath) #TODO dont remove.. do other thing..!
        with open(filepath, 'w') as f:
            f.write(data)

    """
    TODO add doc
    """
    def _refactor( self ):
        print('Running a standard PCA...')
        pca_out1 = pca.PCA(self.meth_data.data.transpose()) 

        print('Compute a low rank approximation of input data and rank sites...')
        x = self._low_rank_approximation(pca_out1.P, pca_out1.U, self.k)
        
        An = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(self.meth_data.data.transpose()).transform(self.meth_data.data.transpose()) #TODO move transpose out?
        Bn = preprocessing.StandardScaler( with_mean = True, with_std = False ).fit(x).transform(x)
        # normalization
        An = An * ( 1 / sqrt((An**2).sum(axis=0)) ) 
        Bn = Bn * ( 1 / sqrt((Bn**2).sum(axis=0)) )

        # find the distance of each site from its low rank approximation.
        distances = self._euclidean_distance(An, Bn)

        ranked_list = distances.argsort()

        print('Compute ReFACTor components...')
        sites = ranked_list[0:self.t]

        pca_out2 = pca.PCA(self.meth_data.data[sites,:].transpose())
        score = pca_out2.P

        print('Saving a ranked list of the data features...')
        data = '\n'.join(['%s\t%s'% (index, self.meth_data.cpgnames[index]) for index in ranked_list])
        self._write_file(self.ranked_output_filename, data)

        print('Saving the ReFACTor components...')
        data = '\n'.join(['\t'.join([str(i) for i in line]) for line in score[:,0:self.k]])
        self._write_file(self.components_output_filename, data)
        
        return score[:,0:self.num_components], ranked_list, pca_out1.P[:,0:self.k]


    def _low_rank_approximation(self, A, B, i):
        """
        calculates low rank approximation between the first i columns of two n X m matrixes A and B 
        Note: Both A and B dimensions are n X m 
        """
        return dot(A[:,0:i], B[:,0:i].transpose())

    def _euclidean_distance(self, A, B):
        """
        calculates euclidean distance between two n X m matrixes A and B
        Note: Both A and B dimensions are n X m 
        """
        return sqrt(((A - B)**2).sum(axis=0))
