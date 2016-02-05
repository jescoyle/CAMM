## This script runs the community assembly of mutualists model (CAMM)

library(dgof) # ks.test for discrete distributions
library(abind) # for growing arrays during simulation run

## Set options, load parameter values and simulation functions
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM'
code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'

# When running on cluster:
working_dir = './Results/'
code_dir = './'

options(stringsAsFactors=F)
source(paste(code_dir,'simulation_functions.R', sep=''))

setwd(working_dir)



