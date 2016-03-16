import os
import sys
import argparse
from refactor_lib import methylation_data
from refactor_lib import refactor_lib


class kEstimateArgumentParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

def run():

        parser = kEstimateArgumentParser(prog=os.path.basename(sys.argv[0]), description = "estimate_k: a tool for estimating the parameter k for ReFACTor by computing a score for the first several eigenvalues of the data matrix. k should be selected to be the number of substantial eigenvalues (in terms of their scores).", epilog = "", add_help=False)
        required = parser.add_argument_group('required arguments')
        required.add_argument('--datafile', required = True, type = str, help = "A data matrix file of beta-normalized methylation levels; see the instructions of refactor.py for more details")
      
        optional = parser.add_argument_group('optional arguments')
        optional.add_argument('-h', '--help', action='help')
        optional.add_argument('--max_k', type = int, default = 10, help = "The maximal eigenvalue to display (DEFAULT=10)")
              
        args = parser.parse_args()

        # load methylation data file
        meth_data = methylation_data.MethylationData(datafile = args.datafile)
        n = meth_data.data.shape[1] # number of samples
        
        # Validate input
        if (args.max_k < 3 or args.max_k > n-1):
            print("max_k must be greater than 2 and smaller than the number of samples")
            sys.exit(2)
        
        refactor_lib.Refactor.estimate_k(methylation_data = meth_data, max_k = args.max_k)
        

if __name__ == "__main__":
  run()
