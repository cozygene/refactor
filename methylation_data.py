import os
import sys
import copy
import logging
from pickle import dump
from numpy import loadtxt, delete, isnan, nanvar, where
from numpy.ma import average, masked_array

COMPRESSED_FILENAME = "methylation_data"

class MethylationData( object ):
    """
    TODO add class doc here
    """
    def __init__(self, datafile):
        data = self._load_and_validate_file_of_dimentions(datafile, 2) 
        self.samples_ids = data[0,:][1:]                 # extract samples ID
        self.cpgnames = data[:,0][1:]                    # extract methylation sites names

        # remove sample ID and sites names from matrix
        # that kind of assignment will create a copy of O[1:,1:], TODO: do we need a new copy here? i don't think so
        # Note that assignment like self.O = O will not create a copy
        self.data = data[1:,1:].astype(float) 
        self.sites_size, self.samples_size = self.data.shape


        logging.debug("Got methylation data with %s sites and %s samples id" % (self.sites_size, self.samples_size))


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
        logging.info("Loading file %s..." % filepath)
        data = loadtxt(filepath, dtype = str)#, converters = lambda x: x if x != 'NA' else 'nan')#,delimiter=';', missing_values='NA', filling_values=nan)# = lambda x: x if x != 'NA' else nan)#, missing_values = '???', filling_values = 0)
        # data = genfromtxt(args.datafile, dtype = str , delimiter=';', usemask = 'True', missing_values = 'NA', filling_values = "???")

        if len(data.shape) != dim:
            logging.error("The file '%s' is not a %sd matrix" % (filepath, dim))
            sys.exit(2)

        return data

    def save(self, methylation_data_filename):
        """
        serializes this object and saves it to methylation_data_filename
        assumes that methylation_data_filename is a valid file 
        """
        with open(methylation_data_filename, 'wb') as f:
            logging.info("Saving methylation data as glint format at %s" % methylation_data_filename)
            dump(self, f)

    def copy(self):
        """
        returns a copy of the object
        """
        return copy.deepcopy(self)
