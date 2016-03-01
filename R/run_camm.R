## This script runs the community assembly of mutualists model (CAMM) on a set of parameter files in its current directory

options(stringsAsFactors=F)
library(reshape) # melt()

# The following arguments should be passed (in order)
# number of cores allocated to this run : MUST BE SUPPLIED!!
# directory with parameter files : defaults to current directory
# directory with simulation scripts : defaults to current directory
# directory in which to save results : defaults to './Results'
args = commandArgs(trailingOnly=T)
print(args)

# Set directories
parm_dir = ifelse(is.na(args[2]), './', args[2])
sim_dir = ifelse(is.na(args[3]), './', args[3])
results_dir = ifelse(is.na(args[4]), './Results/', args[4])

# Load simulation scripts
source(paste0(sim_dir, 'simulation_functions.R'))
source(paste0(sim_dir, 'control_functions.R'))

# Read in parameter files
# Parameter files are designated by starting with 'p_'
file_list = list.files(parm_dir, '^p_')

# Set number of cores
ncores = as.numeric(args[1])

# Define options for simulation
nruns = 100 # 1000 # Number of distinct starts
nchains = 10 # 100 # Number of replicate simulations starting from the same initial metacommunity
sim_mode = 'fixed' # Stopping rule: stop after nreps timesteps
reps = 1000 #1000 #c(1000,2000,4000,8000) # 5000 # Number of timesteps until stop
sim_parms = list(sim_mode=sim_mode, reps=reps)

# Run set of simulations on each parameter
for(f in file_list){
	
	# Read in parameters
	parm_file = paste0(parm_dir, f)
	source(parm_file)

	# Run CAMM
	model_out = run_camm_N(sim_dir, parm_file, nruns, nchains, ncores, sim_parms, simID=runID, save_start=T, save_sim=T, save_dir=results_dir) 
		
	# Analyze results
	S_a = melt(summarize_camm(model_out, 'S', 'a'))
	S_b = melt(summarize_camm(model_out, 'S', 'b'))
	N_comm = melt(summarize_camm(model_out, 'N'))
	comm_summary = rbind(S_a, S_b, N_comm)
	comm_summary = cast(comm_summary, ... ~ response)
	
	cor_a = melt(summarize_camm(model_out, 'cor','a'))
	names(cor_a)[names(cor_a)=='value'] = 'cor_a'
	cor_b = melt(summarize_camm(model_out, 'cor','b'))
	names(cor_b)[names(cor_b)=='value'] = 'cor_b'
	cor_summary = merge(cor_a, cor_b)
	
	write.csv(comm_summary, file=paste0(results_dir, runID, '_comm_summary.csv'), row.names=F)
	write.csv(cor_summary, file=paste0(results_dir, runID, '_cor_summary.csv'), row.names=F)

}

quit(save='no')