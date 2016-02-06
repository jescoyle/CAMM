## This script runs the community assembly of mutualists model (CAMM) on a set of parameter files in its current directory

parm_files = iscrete distributions
library(abind) # for growing arrays during simulation run

## Set options, load parameter values and simulation functions
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM'
code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'

# When running on cluster:
working_dir = './Results/'
code_dir = './'

options(stringsAsFactors=F)
source(paste(code_dir,'simulation_functions.R', sep=''))
source(paste0(code_dir, 'control_functions.R'))

setwd(working_dir)


test1 = run_camm_N('./GitHub/R/', './GitHub/R/parameter_file.R', nruns=4, nchains=10, nparallel=2,
	sim_parms, simID='testone', save_start=T, save_sim=T) 

summarize_camm(test1, 'S', 'a')
summarize_camm(test1, 'S', 'b')
summarize_camm(test1, 'N')
summarize_camm(test1, 'cor', 'a')
summarize_camm(test1, 'cor', 'b')

parm_list = make_parmlist()

write_parms(parm_list, 'testparms', './')