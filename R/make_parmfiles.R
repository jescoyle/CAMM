## This script makes sets of parameter values for CAMM to be run on the cluster

options(stringsAsFactors=F)
setwd('./UNC/Projects/CAMM')

# Set directories
sim_dir = './GitHub/R/'
parm_dir = './Parms/Parms_5/'

dir.create(parm_dir)

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


## RUN 3 ##

# Define parameter sets
breadth_vec = c(0.25, 0.5, 1, 2, 4)
topo_vec = c('one2one','one2many','many2many')
filter_vec = c('opposite','same','none','all')
combos = expand.grid(topo_vec, filter_vec)
combos_sig = expand.grid(breadth_vec, breadth_vec)

for(j in 1:nrow(combos_sig)){
	sigmaA = combos_sig[j,1]
	sigmaB = combos_sig[j,2]

	this_dir = paste0(parm_dir,'sigA_', sigmaA, '-sigB_', sigmaB, '/')
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
			sigma_a1 = 8
			sigma_a2 = sigmaA
			sigma_b1 = sigmaB
			sigma_b2 = 8
		}

		if(envfilt=='same'){
			sigma_a1 = 8
			sigma_a2 = sigmaA
			sigma_b1 = 8
			sigma_b2 = sigmaB
		}

		if(envfilt=='none'){
			sigma_a1 = 8
			sigma_a2 = 8
			sigma_b1 = 8
			sigma_b2 = 8
		}

		if(envfilt=='all'){
			sigma_a1 = sigmaA
			sigma_a2 = sigmaA
			sigma_b1 = sigmaB
			sigma_b2 = sigmaB
		}

		runID = paste0('sigA-', sigmaA, '_sigB-',sigmaB,'_topo-',topology,'_envfilt-',envfilt)

		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}

for(i in breadth_vec){
	windows()
	nicheparms_a$sigma[1] = i
	niches_a = make_niches(S_a, nicheparms_a)
	plot_niches(niches_a, matrix(rep(c(-3,3),2), ncol=2))
}



### RUN 4 : run time ###
topo_vec = c('one2one','one2many','many2many')
filter_vec = c('opposite','same','none','all')
combos = expand.grid(topo_vec, filter_vec)

for(envfilt in filter_vec){
	this_dir = paste0(parm_dir, envfilt, '/')
	dir.create(this_dir)

	for(topology in topo_vec){
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

		runID = paste0('topo-',topology,'_envfilt-',envfilt)

		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}


### RUN 5: species richness ###

topo_vec = c('one2many','many2many')
filter_vec = c('opposite','same')
Sb_vec = 5*2^(0:3)
skew_vec = c(1,2,3,5,10)
combos = expand.grid(Sb_vec, topo_vec)


for(skew in skew_vec){
for(envfilt in filter_vec){
	this_dir = paste0(parm_dir, envfilt, '-', skew, '/')
	dir.create(this_dir)
	
	for(i in 1:nrow(combos)){
		S_b = combos[i,1]
		S_a = skew*S_b
		topology = as.character(combos[i,2])

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
					
		if(topology=='one2many'){
			N_L = S_a
			runID = paste0('topo-',topology,'_envfilt-',envfilt,'_Sa-',S_a,'_Sb-',S_b,'_NL-',N_L)
			parm_list = make_parmlist()
			write_parms(parm_list, paste0('p_', runID), this_dir)
		}

		if(topology=='many2many'){
			links_vec = floor(S_a*c(1.25, 1.5, 2, 3))
			links_vec = links_vec[(links_vec<=S_a*S_b)]

			for(N_L in links_vec){
				runID = paste0('topo-',topology,'_envfilt-',envfilt,'_Sa-',S_a,'_Sb-',S_b,'_NL-',N_L)
				parm_list = make_parmlist()
				write_parms(parm_list, paste0('p_', runID), this_dir)
			}
		}
	}
}}







