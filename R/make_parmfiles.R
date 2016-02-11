## This script makes sets of parameter values for CAMM to be run on the cluster

options(stringsAsFactors=F)
setwd('./UNC/Projects/CAMM')

# Set directories
sim_dir = './GitHub/R/'
parm_dir = './Parms/'


# Read in base parameter file
source(paste0(sim_dir, 'parameter_file.R'))

# Read in functions
source(paste0(sim_dir, 'control_functions.R'))


## RUN 1 ##
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


## RUN 2 ##

# Define parameter sets
omega_vec = c(0, .5, .9, 1)
topo_vec = c('one2one','one2many','many2many')
filter_vec = c('opposite','same','none', 'all')
combos = expand.grid(topo_vec, filter_vec)

for(o in omega_vec){
	omega = o
	this_dir = paste0(parm_dir,'omega_', o, '/')
	dir.create(this_dir)

	for(i in 1:nrow(combos)){
		topology = as.character(combos[i,1])
		
		if(topology=='one2one'){
			S_a = 30
			S_b = 30
			N_L = 30
		}
		
		if(topology=='one2many'){
			S_a = 30
			S_b = 10
			N_L = 30
		}

		if(topology=='many2many'){
			S_a = 30
			S_b = 10
			N_L = 60
		}

		envfilt = combos[i,2]
	
		if(envfilt=='opposite'){
			sigma_a1 = 10
			sigma_a2 = 0.5
			sigma_b1 = 0.5
			sigma_b2 = 10
		}

		if(envfilt=='same'){
			sigma_a1 = 10
			sigma_a2 = 0.5
			sigma_b1 = 10
			sigma_b2 = 0.5
		}

		if(envfilt=='none'){
			sigma_a1 = 10
			sigma_a2 = 10
			sigma_b1 = 10
			sigma_b2 = 10
		}

		if(envfilt=='all'){
			sigma_a1 = 0.5
			sigma_a2 = 0.5
			sigma_b1 = 0.5
			sigma_b2 = 0.5
		}

		runID = paste0('o-',o,'_topo-',topology,'_envfilt-',envfilt)

		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}














