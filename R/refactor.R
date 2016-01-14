
if(!exists("refactor", mode="function")) source("refactor_lib.R")# args <- commandArgs(trailingOnly = TRUE)

args <- commandArgs(trailingOnly = TRUE)
data_file=args[1] 
k=as.numeric(args[2]) 
num_comp=args[3] 

# run refactor
output <- refactor(data_file, k, num_components=num_comp)
