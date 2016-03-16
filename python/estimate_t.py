import os
import sys
import argparse
from refactor_lib import methylation_data
from refactor_lib import refactor


class tEstimateArgumentParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)

def run():
        parser = tEstimateArgumentParser(prog=os.path.basename(sys.argv[0]), description = "estimate_t: a tool for estimating the parameter t for ReFACTor by computing a score for the sites based on their low-rank approximation. t should be selected to be the number of sites with substantial signal.", epilog = "", add_help=False)
        required = parser.add_argument_group('required arguments')
        required.add_argument('--datafile', required = True, type = str, help = "A data matrix file of beta-normalized methylation levels; see the instructions of refactor.py for more details")
        required.add_argument('--k', required = True, type = int, help = "The number of assumed cell types")
      
        optional = parser.add_argument_group('optional arguments')
        optional.add_argument('-h', '--help', action='help')
        optional.add_argument('--numsites', type = int, default = 5000, help = "The number of sites to display (DEFAULT=min{5000,#sites})")
              
        args = parser.parse_args()

        # load methylation data file
        meth_data = methylation_data.MethylationData(datafile = args.datafile)

        m = meth_data.data.shape[0] # number of sites
        if (m < args.numsites):
            args.numsites = m

        # Validate input
        if (args.numsites < 1):
            print("numsites must be greater than 0")
            sys.exit(2)

        refactor.Refactor.estimate_t(methylation_data = meth_data, k = args.k, numsites = args.numsites)
        

if __name__ == "__main__":
  run()
    

