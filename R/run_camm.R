## This script runs the community assembly of mutualists model (CAMM)

library(abind)

## Set options, load parameter values and simulation functions
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM'
code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'
setwd(working_dir)

options(stringsAsFactors=F)
source(paste(code_dir,'parameter_file.R', sep=''))
source(paste(code_dir,'simulation_functions.R', sep=''))

## Instantiate communities and mutualistic network
topo = make_topo(S_a, S_b, N_L, topology) # S_a x S_b matrix of association probabilities
topo_names = name_topo(topo) # S_a x S_b matrix with integers labeling specific associations
sites = make_sites(N_C, rho_z) # N_C x 2 matrix of env values
comm = make_comm(N_C, N, N_L) # initial N_C x N matrix of integers indicating which association is present
#comm = rand_comm(N_C, N, N_L) # Make random community for testing simulation
poolA = calc_pool(comm, topo_names, 1) # N_C x S_a matrix of presence of mutualist a in communities based on associations present
poolB = calc_pool(comm, topo_names, 2) # N_C x S_b matrix of presence of mutualist b in communities based on associations present

## Generate random niches for each mutualist
niches_a = make_niches(S_a, nicheparms_a) # array of S_a x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
niches_b = make_niches(S_b, nicheparms_b)	# array of S_b x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)

## Generate global species abundance distribution
# Define niche quality with which global abudance should be correlated
if(gsad_cond_a$type=='none'){
	gsad_a = make_gsad(S_a, gsad_dist_a)
} else {
	if(gsad_cond_a$type=='breadth'){
		cond_a = apply(niches_a[,'sigma',], 1, prod)
	}
	if(gsad_cond_a$type=='links'){
		cond_a = rowSums(topo)
	}
	gsad_a = make_gsad(S_a, gsad_dist_a, cond_a, gsad_cond_a$rho)
}

if(gsad_cond_b$type=='none'){
	gsad_b = make_gsad(S_b, gsad_dist_b)
} else {
	if(gsad_cond_b$type=='breadth'){
		cond_b = apply(niches_b[,'sigma',], 1, prod)
	}
	if(gsad_cond_b$type=='links'){
		cond_b = colSums(topo)
	}
	gsad_b = make_gsad(S_b, gsad_dist_b, cond_b, gsad_cond_b$rho)
}

# Plot niches for visualization
# plot_niches(niche_array, matrix of x-lims, site_env)
#plot_niches(niches_a, matrix(c(-3,3,-3,3), 2, 2), add_env=sites)
#plot_niches(niches_b, matrix(c(-3,3,-3,3), 2, 2), add_env=sites)

# Plot topology for visualization
#plot_topo(topo, orderby='degree')

# Calculate N_C x N_L matrix of association probabilities at each site
assoc_probs = matrix(0, N_C, N_L)
for(i in 1:N_C){
for(j in 1:N_L){
	env = sites[i,]

	partners = which(topo_names==j, arr.ind=T)
	n_a = niches_a[partners[1],'mu',]
	n_b = niches_b[partners[2],'mu',]

	assoc_probs[i,j] = assoc_prob_func(env, n_a, n_b)
}}


## A simulation of Markov process: dispersal implicit, rates independent for both mutualists


### Run simulation for a given number of time steps
if(sim_mode=='fixed'){

# Generate empty array to hold changes in community and add initial community
comm_records = array(NA, c(N_C, N, reps+1)) 
comm_records[,,1] = comm
poolA_records = array(NA, c(N_C, S_a, reps+1)) 
poolA_records[,,1] = poolA
poolB_records = array(NA, c(N_C, S_b, reps+1))
poolB_records[,,1] = poolB
Tmat_records = array(NA, c(N_C, N_L+1, N_L+1, reps+1))

## Run simulation
for(step in 1:reps){
	# Mutualist mortality: mutualists present as associations cannot die
	poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
	poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)

	# Mutualists disperse into communities independently and with probability based on niche-based filtering	
	poolA = disperse(sites, niches_a, poolA, gsad_a)
	poolB = disperse(sites, niches_b, poolB, gsad_b)

	# Calculate transition matrices for each site
	T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs) # array of square matrices of transition probabilities from one association to another
	Tmat_records[,,,step] = T_mat

	# Identify current state of each space and transition based on random numbers
	new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
	for(i in 1:nrow(comm)){
	for(j in 1:ncol(comm)){
		new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
	}}
	comm = new_comm

	# Save communities
	comm_records[,,step+1] = comm
	poolA_records[,,step+1] = poolA
	poolB_records[,,step+1] = poolB

}
} # CLOSES if(sim_mode=='fixed') 



### Run simulation until N chains converge on equilibrial community dynamics
if(sim_mode=='converge'){

# WHAT CRITERIA COUNTS AS CONVERGENCE?
# CONVERGENCE = CHAINS ARE SAMPLING FROM THE SAME DISTRIBUTION OF COMMUNITIES 
# from Brooks and Gelman 1998 : CONVERGENCE = WITHIN-CHAIN PERCENTILE = MIXTURE OF SEQUENCES PERCENTILE, AND BOTH STABILIZE


# MAY WANT TO RE-WRITE THIS SECTION TO RUN IN PARALLEL

# Split community into nchain identical instances
# COULD CHANGE THIS TO START FROM RANDOM COMMUNITIES
comm_arr = array(comm, dim=c(nrow(comm), ncol(comm), nchains))
poolA_arr = array(poolA, dim=c(nrow(poolA), ncol(poolA), nchains))
poolB_arr = array(poolB, dim=c(nrow(poolB), ncol(poolB), nchains))

# Generate empty array to hold changes in community and add initial community
comm_records = array(NA, c(N_C, N, reps+1, nchains)) 
comm_records[,,1,] = comm_arr
poolA_records = array(NA, c(N_C, S_a, reps+1, nchains)) 
poolA_records[,,1,] = poolA_arr
poolB_records = array(NA, c(N_C, S_b, reps+1, nchains))
poolB_records[,,1,] = poolB_arr
Tmat_records = array(NA, c(N_C, N_L+1, N_L+1, reps+1, nchains))

## Run simulation until convergence
window_size = 50 # Number of iterations to run before calculating convergence criteria
burnin = 100
step = 1

# Make vector to hold convergence flags

while((CONVERGENCE CRITERION NOT MET)&(step > window_size)){
	for(step in 101:200){
	for(k in 1:nchains){
		comm = comm_arr[,,k]
		poolA = poolA_arr[,,k]
		poolB = poolB_arr[,,k]

		# Mutualist mortality: mutualists present as associations cannot die
		poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
		poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)

		# Mutualists disperse into communities independently and with probability based on niche-based filtering	
		poolA = disperse(sites, niches_a, poolA, gsad_a)
		poolB = disperse(sites, niches_b, poolB, gsad_b)

		# Calculate transition matrices for each site
		T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs) # array of square matrices of transition probabilities from one association to another
		Tmat_records[,,,step, k] = T_mat

		# Identify current state of each space and transition based on random numbers
		new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
		for(i in 1:nrow(comm)){
		for(j in 1:ncol(comm)){
			new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
		}}
		comm = new_comm

		# Save communities
		comm_arr[,,k] = comm
		poolA_arr[,,k] = poolA
		poolB_arr[,,k] = poolB
		comm_records[,,step+1,k] = comm
		poolA_records[,,step+1,k] = poolA
		poolB_records[,,step+1,k] = poolB
	}
	}
	# Calculate convergence criteria
	if((step %% window_size == 0)&(step > burnin)){
		# Determine observations to use in calculation of R_hat
		n_obs = floor(length((burnin+1):step)/2)
		use_obs = (step-n_obs+1) : step

		# Calculate R_hat for species richness
		commstats = array(NA, dim=c(n_obs, nchains, N_C, 4), 
			dimnames=list(step=use_obs, chain=1:nchains, comm = 1:N_C, statistic=c('S_species','S_a','S_b','N')))
		for(n in use_obs){
		for(k in 1:dim(comm_records)[4]){	
			commstats[as.character(n),k,,] = as.matrix(calc_commstats(comm_records[,,n,k], topo_names))
		}}	

		

		

		# COMPARE COMPOSITION~ENV CORRELATIONS?

	}
	step = step + 1

} # CLOSES while(converegence) loop

} # CLOSES if(sim_mode=='converge')

