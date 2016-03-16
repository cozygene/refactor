#!/usr/bin/env python

import os
import sys
import argparse
from refactor_lib import refactor

REFACTOR_FORMATTED_EXTENSION = ".refactor"

class RefactorArgumentParser(argparse.ArgumentParser):

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)


def run():
      parser = RefactorArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description = "ReFACTor: Reference-Free Adjustment for Cell-Type composition",
                                     epilog = "",
                                      add_help=False)
      required = parser.add_argument_group('required arguments')
      required.add_argument('--k',        required = True, type = int,                help = "The number of assumed cell types")
      required.add_argument('--datafile', required = True, type = str,                help = "A data matrix file of beta-normalized methylation levels; see the README file for more details")
      
      optional = parser.add_argument_group('optional arguments')
      optional.add_argument('-h', '--help', action='help')
      optional.add_argument('--covarfile',                       type = str,   help = "A covariates file; see the README file for more details")
      optional.add_argument('--t',                         type = int, default = 500, help = "The number of sites to use for computing the ReFACTor components (DEFAULT=500)")
      optional.add_argument('--numcomp',                   type = int,                help = "The number of ReFACTor components to output (DEFAULT=k)")
      optional.add_argument('--stdth',                   type = float, default = 0.02, help = "The standard deviation threshold for excluding low variance sites (DEFAULT=0.02)")
      optional.add_argument('--out',                       type = str, default ="refactor",   help = "Changes the prefix of the output files (DEFAULT='refactor')")
      
      args = parser.parse_args()

      return refactor.Refactor(datafile = args.datafile, 
                                k = args.k, 
                                covarfile = args.covarfile,
                                t = args.t, 
                                num_components = args.numcomp,
                                stdth = args.stdth,
                                out = args.out)


if __name__ == "__main__":
  run()
    
