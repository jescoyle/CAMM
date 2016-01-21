## This script runs the community assembly of mutualists model (CAMM)

library(dgof) # ks.test for discrete distributions
library(abind) # for growing arrays during simulation run

## Set options, load parameter values and simulation functions
#working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM'
#code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'

# When running on cluster:
working_dir = './Results/'
code_dir = './'

options(stringsAsFactors=F)
source(paste(code_dir,'parameter_file.R', sep=''))
source(paste(code_dir,'simulation_functions.R', sep=''))

setwd(working_dir)

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

# Save community parameters
save(topo, topo_names, sites, niches_a, niches_b, gsad_a, gsad_b, assoc_probs, file=paste('sim_object_', runID, '.RData', sep=''))


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
Tmat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs)
Tmat_arr = array(Tmat, dim=c(dim(Tmat), nchains))
Rhat_arr = array()

# Generate empty array to hold changes in community and add initial community
comm_records = array(comm_arr, dim=c(dim(comm_arr), 1)) 
poolA_records = array(poolA_arr, dim=c(dim(poolA_arr),1)) 
poolB_records = array(poolB_arr, dim=c(dim(poolB_arr),1))
Tmat_records = array(Tmat_arr, dim=c(dim(Tmat_arr),1))


## Run simulation until convergence
window_size = 50 # Number of iterations to run before re-calculating convergence criteria
burnin = 100
thin = 1
Rhat_tol = 0.1 # Tolerance for R_hat statistic to be considered converged

# Set initial counters
converged = F
step = 1

# COMMUNITIES ARE NOT CONVERGING UNDER THESE CRITERIA- MAY WANT TO RECONSIDER
while(!converged){
	print(paste('Step =', step))

	# Do one step through the simulation
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
		Tmat_arr[,,,k] = T_mat

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
	}
	
	# Save communities from this step
	Tmat_records = abind(Tmat_records, Tmat_arr)
	comm_records = abind(comm_records, comm_arr)
	poolA_records = abind(poolA_records, poolA_arr)
	poolB_records = abind(poolB_records, poolB_arr)
	
	# Calculate convergence criteria if past burnin and at correct re-test interval
	if((step %% window_size == 0)&(step > burnin)){
		# Determine observations to use in calculation of R_hat
		window = floor(length((burnin+1):step)/2)
		use_obs = seq((step-window+1), step, thin)
		n_obs = length(use_obs)
		
		print(paste('Number of obs =', window, '/', thin))
		
		# RICHNESS DOESN'T APPEAR TO STABILIZE ENOUGH FOR K-S TEST TO CONCLUDE CONVERGENCE
		# Calculate community statistics for each chain and each community
		#commstats = array(NA, dim=c(n_obs, nchains, N_C, 4), 
		#	dimnames=list(step=use_obs, chain=1:nchains, comm = 1:N_C, statistic=c('S_species','S_a','S_b','N')))
		#for(n in use_obs){
		#for(k in 1:nchains){	
		#	commstats[as.character(n),k,,] = as.matrix(calc_commstats(comm_records[,,k,n], topo_names))
		#}}	
		
		# For each community statistic, use K-S test to determine whether each chain is sampling from empirical distribution function of all chains mixed
		#rich_flags = sapply(dimnames(commstats)$statistic, function(y){
		#	sapply(1:N_C, function(i){
		#		pvals = sapply(1:nchains, function(k) ks.test(commstats[,k,i,y],ecdf(commstats[,,i,y]))$p.value)
		#		prod(pvals > 0.95)				
		#	})
		#}) > 0
		
		#print('Community Statistics Convergence:')
		#print(colSums(rich_flags)/nrow(rich_flags))
	
		# Calculate community composition~env correlation for each chain
		corrstats = array(NA, dim=c(n_obs, nchains, 3, ncol(sites)), 
			dimnames=list(step=use_obs, chain=1:nchains, community=c('species','a','b'), env=1:ncol(sites)))
		for(n in use_obs){
		for(k in 1:nchains){
			corrstats[as.character(n), k, , ] = calc_envcorr(comm_records[,,k,n], topo_names, sites, c(), binary=F)
		}}

		# Calculate R_hat among chains
		Rhats = sapply(dimnames(corrstats)$community, function(y){
			sapply(1:ncol(sites), function(x){
				calc_Rhat(corrstats[,,y,x], alpha=.95)
			})
		})

		# Save records of Rhat
		if(exists('Rhat_records')){ 
			Rhat_records = abind(Rhat_records, Rhats)
		} else {
			Rhat_records = array(Rhats, dim=c(dim(Rhats),1), 
				dimnames=list(env=1:dim(Rhats)[1], community=dimnames(Rhats)[[2]]))
		}

		corr_flags = abs(Rhats-1) < Rhat_tol
		
		print('Community~Env Correlation Convergence:')
		print(corr_flags)
	
		# Test whether all flags are TRUE
		#converged = prod(c(corr_flags, rich_flags))>0
		converged = prod(corr_flags)>0

		# Save results so far
		save(comm_records, poolA_records, poolB_records, Tmat_records, Rhat_records, 
			file=paste('sim_results_',runID,'.RData',sep=''))

		# WORKING HERE
		# Plot convergence criteria
		# LOOK INTO BINOMIAL DISSIMILARITY IN VEGDIST FOR COMPARING COMMUNITIES ACROSS CHAINS

	}
	step = step + 1

} # CLOSES while(converegence) loop

} # CLOSES if(sim_mode=='converge')

# Close R and don't save the session
quit('no')

