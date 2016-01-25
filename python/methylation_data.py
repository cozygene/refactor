import os
import sys
from numpy import loadtxt, delete, isnan, nanvar, where, std
from numpy.ma import average, masked_array
import copy

COMPRESSED_FILENAME = "methylation_data"

class MethylationData( object ):

    def __init__(self, datafile):
        data = self._load_and_validate_file_of_dimentions(datafile, 2) 
        self.samples_ids = data[0,:][1:]                 # extract samples ID
        self.cpgnames = data[:,0][1:]                    # extract methylation sites names

        self.data = data[1:,1:].astype(float) 
        self.sites_size, self.samples_size = self.data.shape


    def _validate_file_path(self, filepath):
        if not os.path.exists(filepath) :
            print("ERROR: The file '%s' doesn't exist. Exiting" % filepath)
            sys.exit(2)

    def _load_and_validate_file_of_dimentions(self, filepath, dim):
        """
        validates that a file exists and that it is a matrix from dimentions dim.
        """
        if filepath is None:
            return None

        self._validate_file_path(filepath)
        print("Loading file %s..." % filepath)
        data = loadtxt(filepath, dtype = str)

        if len(data.shape) != dim:
            print("ERROR: The file '%s' is not a %sd matrix" % (filepath, dim))
            sys.exit(2)

        return data

    def _filter_sites_by_std(self, th):
        """
        Removes sites with std lower than the specified threshold.
        """
        stds = std(self.data,1)
        stds_sorted = sorted(stds)
        stds_sorted_ind = stds.argsort()
        p = self.sites_size-1
        while (p > 0):
            if stds_sorted[p] < th:
                p += 1
                break
            p -= 1
        if (p == self.sites_size):
            print("ERROR: the provided stdth parameter excludes all sites")
            sys.exit(2)
        sites = stds_sorted_ind[p:]
        self.data = self.data[sites,:]
        # Update class fields
        self.sites_size, self.samples_size = self.data.shape
        self.cpgnames = self.cpgnames[sites]

    def _copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)

