import os
import sys
import argparse
from refactor_lib import Refactor
from methylation_data import MethylationData

REFACTOR_FORMATTED_EXTENSION = ".refactor"

class RefactorArgumentParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)


def run():
      parser = RefactorArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description = "ReFACTor: Reference-Free Adjustment for Cell-Type composition",
                                     epilog = "")
      
      parser.add_argument('--datafile', required = True, type = str,                help = "A data matrix file of beta-normalized methylation levels; see the README file for more details")
      parser.add_argument('--k',        required = True, type = int,                help = "The number of assumed cell types")
      parser.add_argument('--t',                         type = int, default = 500, help = "The number of sites to use for computing the ReFACTor components (DEFAULT=500)")
      parser.add_argument('--numcomp',                   type = int,                help = "The number of ReFACTor components to output (DEFAULT=K)")
      parser.add_argument('--out',                       type = str, default ="",   help = "changes the prefix of the output file ")

      args = parser.parse_args()
      
      # load methylation data file
      meth_data = MethylationData(datafile = args.datafile)
      
      return Refactor(methylation_data = meth_data, 
                                k = args.k, 
                                t = args.t, 
                                num_components = args.numcomp, 
                                ranked_output_filename = args.out + Refactor.RANKED_FILENAME, 
                                components_output_filename  = args.out + Refactor.COMPONENTS_FILENAME)


if __name__ == "__main__":
  run()
    
