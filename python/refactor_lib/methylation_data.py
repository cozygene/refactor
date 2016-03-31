import os
import sys
from numpy import loadtxt, delete, isnan, nanvar, where, std
from numpy.ma import average, masked_array
import copy
from bisect import bisect_right

COMPRESSED_FILENAME = "methylation_data"

class MethylationData( object ):

    def __init__(self, datafile):
        data = self._load_and_validate_file_of_dimensions(datafile, 2) 
        self.samples_ids = data[0,:][1:]                 # extract samples ID
        self.cpgnames = data[:,0][1:]                    # extract methylation sites names

        self.data = data[1:,1:].astype(float) 
        self.sites_size, self.samples_size = self.data.shape


    def _validate_file_path(self, filepath):
        if not os.path.exists(filepath) :
            print("ERROR: The file '%s' doesn't exist. Exiting" % filepath)
            sys.exit(2)

    def _load_and_validate_file_of_dimensions(self, filepath, dim):
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

    def _exclude_sites_from_data(self, sites_indicies_list):
        """
        this function removes from the data the cpg sites which indices found in sites_indicies_list
        it updates the sites_size, the cpgnames list and the list holds the average value per site
        """
        self.data = delete(self.data, sites_indicies_list, axis = 0)
        self.cpgnames = delete(self.cpgnames, sites_indicies_list)
        self.sites_size = len(self.cpgnames)

    def _filter_sites_by_std(self, th):
        """
        Removes sites with std lower than the specified threshold.
        """
        stds = std(self.data,1)
        stds_sorted_ind = stds.argsort()
        stds_sorted = stds[stds_sorted_ind]

        p = bisect_right(stds_sorted, th) # binary search: find leftmost value greater than th and returns it's index
        if (p == self.sites_size):
            print("ERROR: the provided stdth parameter excludes all sites")
            sys.exit(2)
        exclude_ind = stds_sorted_ind[:p]
        self._exclude_sites_from_data(exclude_ind)
      

    def _copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)

