## This script makes sets of parameter values for CAMM to be run on the cluster

options(stringsAsFactors=F)
setwd('./UNC/Projects/CAMM')

# Set directories
sim_dir = './GitHub/R/'
parm_dir = './Parms/'


# Read in base parameter file
source(paste0(sim_dir, 'parameter_file.R'))

# Define parameter sets
omega_vec = seq(0, 1, .1)
mort_rate_a_vec = c(1, 10, 20)
mort_rate_b_vec = c(1, 10, 20)
combos = expand.grid(mort_rate_a_vec, mort_rate_b_vec)

# For each omega, make a new directory of parameter files
for(o in omega_vec){
	omega = o

	this_dir = paste0(parm_dir, 'omega_', o, '/')
	dir.create(this_dir)
	
	for(i in 1:nrow(combos)){
		mort_rate_a = combos[i,1]
		mort_rate_b = combos[i,2]
		runID = paste0('o-',omega,'_mra-',mort_rate_a,'_mrb-',mort_rate_b)  

		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}



