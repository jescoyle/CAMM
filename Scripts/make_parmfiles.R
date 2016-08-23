## This script makes sets of parameter values for CAMM to be run on the cluster

options(stringsAsFactors=F)
setwd('C:/Users/jrcoyle/Documents/Research/CAMM/')

# Load package
library(CAMM)

# Set directories
sim_dir = './GitHub/Scripts/'


dir.create(parm_dir)

# Read in base parameter file
source(file.path(sim_dir, 'parameter_file.R'))

# Set baseline parameters that differ from example parameter file
S_a = 30
S_b = 10
N_C = 20
N = 100

################################################################
### Simulation Experiments



## RUN 1 ##
parm_dir = './Parms/Parms_1'

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
parm_dir = './Parms/Parms_2'

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


## RUN 3 & 3B ##
parm_dir = './Parms/Parms_3'
parm_dir = './Parms/Parms_3B'

# For 3B
omega = 0

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
parm_dir = './Parms/Parms_4'

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
parm_dir = './Parms/Parms_5'

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


### RUN 6: uneven species abundances that are correlated with env niches ###
parm_dir = './Parms/Parms_6'

maxN_vec = 2^(1:5)
r_vec = c(0.5, 0.9)
envfilt_vec = c('opposite','same','none','all')
a_corr = 0:2
b_corr = 0:2

combos_1 = expand.grid(maxN_vec, r_vec)
combos_2 = expand.grid(envfilt_vec, a_corr, b_corr)
combos_2 = combos_2[-(1:4),]

for(i in 1:nrow(combos_1)){
	maxN = combos_1[i,1]	
	r = combos_1[i,2]
	
	this_dir = paste0(parm_dir, maxN, '-', r, '/')
	dir.create(this_dir)

	for(j in 1:nrow(combos_2)){
		source(paste0(sim_dir, 'parameter_file.R'))

		envfilt = as.character(combos_2[j, 1])
		corrA = combos_2[j,2]
		corrB = combos_2[j,3]	
		
		gsad_dist_a$maxN = maxN
		gsad_dist_b$maxN = maxN
		gsad_dist_a$corr[corrA,1] = r
		gsad_dist_b$corr[corrB,1] = r
		
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

		runID = paste0('maxN-', maxN, '_r-', r, '_corrA-', corrA, '_corrB-', corrB, '_envfilt-',envfilt)
		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}

# Make parm file with uncorrelated gsad and niche optima
for(maxN in maxN_vec){
	this_dir = paste0(parm_dir, maxN, '-0/')
	dir.create(this_dir)
	
	gsad_dist_a$maxN = maxN
	gsad_dist_b$maxN = maxN

	for(envfilt in envfilt_vec){
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
	
		runID = paste0('maxN-', maxN, '_r-0_corrA-0_corrB-0_envfilt-',envfilt)
		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}


### RUN 7 ###
parm_dir = './Parms/Parms_7'

## Effect of changing community size

topo_vec = c('one2many','many2many')
envfilt_vec = c('same','opposite')
N_vec = S_a*2^(-1:4)
NC_vec = c(10, 20, 40, 80)
combos = expand.grid(NC_vec, N_vec)
S_a = 30
S_b = 10

for(envfilt in envfilt_vec){
for(topology in topo_vec){
	this_dir = paste0(parm_dir, topology, '-', envfilt, '/')
	dir.create(this_dir)
	
	for(i in 1:nrow(combos)){
		N = combos[i,2]
		N_C = combos[i,1]

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
			N_L = 30
		}

		if(topology=='many2many'){
			N_L = 60
		}

		runID = paste0('topo-',topology,'_envfilt-',envfilt,'_N-',N,'_NC-',N_C)
		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
	}
}}


##################################################################
### Empirical Runs

parm_dir = './Parms/MP'
setwd(parm_dir)


### Myco-Photo TRFLP ###

# Set mycobiont richness
S_a = 57

# Set number of samples
N_C = 54

# Define community size
# Assumes 9 x 15 cm area, where an individual takes up 4 cm2
N = 34

# Define number of photobionts
# Based on observed number of taxa in TRFLP community data matrices under different abdance criteria
# see 'Analysis/Derived_Data/compare_strain_richness.csv'
Sb_vec = c(15, 20, 25, 30, 35, 40, 45, 60, 70, 75, 80, 85, 90, 95, 105)

# Define topology
topology = 'many2many'

# Define factor for determining number of links
NL_vec = c(57, 86, 114) # x 1, 1.5, 2

# Define environmental variables
data_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/Mycobiont - Photobiont/Analysis/Data/'
samples = read.csv(file.path(data_dir, 'samples.csv'))
env = read.csv(file.path(data_dir, 'loggerdata.csv'))
samples = merge(samples, env)
use_env = samples[samples$SampID %in% 19:72, c('Light_mean','Vpd_daysatfreq')] # based on response of communities and correlation structure of env data
std_env = scale(use_env, center=T, scale=T)

# Define filtering
# Hyp: mycos filter on VPD (env2) while photos filter on light (env1)
filt_vec = c('none', 'myco', 'photo', 'light', 'vpd', 'opposite', 'all')

# Define filtering strength
sig = 0.5



# Define niche parameters to match observed data
# Assumes that observed env gradient is small relative to species potential range in niche optima
mu_a1 = 3
mu_a2 = 3
mu_b1 = 3
mu_b2 = 3

prcomp(use_env, center=T, scale=T)

# Make parameter files

for(NLfact in NL_vec){
	
	N_L = S_a*NLfact
	
	for(envfilt in envfilt_vec){

		mkdir(paste0('NL-',S_b))

	

		if(envfilt=='photo'){
			sigma_a1 = 10
			sigma_a2 = 10
			sigma_b1 = sig
			sigma_b2 = 10
		}
		
		if(envfilt=='myco'){
			sigma_a1 = 10
			sigma_a2 = sig
			sigma_b1 = 10
			sigma_b2 = 10
		}

		if(envfilt=='opposite'){
			sigma_a1 = 10
			sigma_a2 = sig
			sigma_b1 = sig
			sigma_b2 = 10
		}

		if(envfilt=='light'){
			sigma_a1 = 0.5
			sigma_a2 = 10
			sigma_b1 = 0.5
			sigma_b2 = 10
		}

		if(envfilt=='vpd'){
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

		
			
		for(S_b in Sb_vec){
			

					runID = paste0('topo-',topology,'_envfilt-',envfilt,'_N-',N,'_NC-',N_C)
		parm_list = make_parmlist()
		write_parms(parm_list, paste0('p_', runID), this_dir)
			



		}
	









}











