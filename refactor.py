"""""
ReFactor
"""""
import os
import sys
import argparse
from refactor_lib import Refactor
from methylation_data import MethylationData
import configurelogging
configurelogging.configureLogging('')
import logging

REFACTOR_FORMATTED_EXTENSION = ".refactor"

class RefactorArgumentParser(argparse.ArgumentParser):#,argparse._ActionsContainer, GlintArgumentGroup):

    # TODO add epilog

    def error(self, message):
        sys.stderr.write('Error: %s\n' % message)
        print("To see the full help: %s -h/--help" % self.prog)
        sys.exit(2)


def run():
      parser = RefactorArgumentParser(prog=os.path.basename(sys.argv[0]),
                                     description = "<< add help before >>",
                                     epilog = "<< add help after >>"))
      
      parser.add_argument('--datafile', required = True, type = str,                help = "A data matrix file of beta-normalized methylation levels or a .glint file")
      parser.add_argument('-k',         required = True, type = int,                help = "The number of assumed cell types")
      parser.add_argument('-t',                          type = int, default = 500, help = "The number of sites to use for computing the ReFACtor components (DEFAULT=500)")
      parser.add_argument('--numcomp',                   type = int,                help = "The number of ReFACTor components to output (DEFAULT=K)")
      parser.add_argument('--out',                       type = str, default ="",   help = "changes the prefix of the output file ")
      parser.add_argument('--gsave',                     action='store_true',       help = "Save the data in a glint format; makes following executions faster")

      args = parser.parse_args()
      
      # load methylation data file
      if args.datafile.endswith(REFACTOR_FORMATTED_EXTENSION):
          with open(args.datafile,'rb') as f:
              logging.info("Loading glint file: %s..." % args.datafile)
              meth_data = load(f)
              logging.debug("Got methylation data with %s sites and %s samples id" % (self.module.sites_size, self.module.samples_size))
      else:
          meth_data = MethylationData(datafile = args.datafile)

      # save serialized methylation data if asked
      if args.gsave:
          meth_data.save(output_perfix + methylation_data.COMPRESSED_FILENAME + REFACTOR_FORMATTED_EXTENSION)

      
      return Refactor(methylation_data = meth_data, 
                                K = args.k, 
                                t = args.t, 
                                num_components = args.numcomp, 
                                ranked_output_filename = args.out + refactor.RANKED_FILENAME, 
                                components_output_filename  = args.out + refactor.COMPONENTS_FILENAME)


if __name__ == "__main__":
  run()
    
