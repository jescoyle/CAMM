#' Run a simulation
#'
#' Runs a simulation on an existing metacommunity
#'
#' Simulations are run for a fixed number of timesteps (\code{reps}) on the 
#' metacommunity provided in \code{metacomm}, or if not provided, on the 
#' metacommunity objects that exist in the environment where the function is 
#' called. The simulation keeps track of the associations (\code{comm}),
#' host species (\code{poolA}), and symbiont species (\code{poolB}) present
#' at each site at each timestep to be saved (\code{save_steps}, defaults to all.
#' It also records the transition matrix used to generate the status of the 
#' at each of these timesteps (\code{Tmat}. These objects are returned as a 
#' list.
#' 
#' For each timestep, the simulation executes the following operations in 
#' order:
#' \enumerate{
#'		\item \strong{Local Extinction} : Stochastic extinction of species from 
#'			sites where they are not established in an association at a microsite. 
#'			See \code{\link{die}}.
#'		\item \strong{Immigration} : Stochastic immigration of species from mainland
#'			into species pools at local sites. The probability that a species becomes
#'			established in a site's local pool is a function of the number of propagules
#'			arriving and the match between the species' niche and environmental
#'			conditions. See \code{\link{disperse}}.
#'		\item \strong{Transition Communities} : Calculate a transition matrix
#'			for each site based on species present in the species pool and then
#'			stochastically change microsite states based on these transition
#'			probabilities. This is the step where established associations can
#'			die and empty microsites are colonized from species present in local
#'			site species pools. See \code{\link{calc_probs}} and 
#'			\code{\link{transition}}.
#'	}
#'
#' @note A previous version of the function allowed the simulation to run until
#' multiple 'chains' reached a convergence criteria. This mode has been removed,
#' but the code remains commented out in the function in case a future user wants
#' to utilize it.
#'
#' @param metacomm list containing objects describing the metacommunity, 
#' 	as returned by \code{\link{initialize_camm}}. If not supplied, then then
#' 	environment must contain the following:
#' 	\describe{
#' 		\item{\code{comm}}{metacommunity object (see \code{make_comm)})}
#' 		\item{\code{poolA} and \code{poolB}}{species pools for hosts and
#' 			symbionts (see \code{calc_pool})}
#' 		\item{\code{niches_a}, \code{niches_b}, \code{gsad_a}, \code{gsad_b}}{
#' 			host ans symbiont niche parameters and global species abundance 
#' 			distributions (see \code{make_niches_gsad})}
#' 		\item{\code{assoc_probs}}{site x association matrix giving the
#' 			probability of each association forming at each site}
#' 		\item{\code{topo} and \code{topo_names}}{binary host x symbiont species
#' 			matrix giving the association network (see \code{make_topo}) and 
#' 			a versionwith integers identifying each association (see 
#' 			\code{name_topo})}
#' 	}
#' @param reps (required) number of timesteps to run simulation
#' @param save_steps vector indicating which timesteps to save. Defaults to all.
#' @return a list containing objects describing states of the metacommunity 
#' 	through time (see details)
#'
#' @export

run_camm = function(metacomm=NULL, reps=NA, save_steps=NA){
	# Previous versions included sim_mode as a parameter (see code commented out below)
	# However, simulations were very slow to converge so the option has been removed and all
	# simulations are performed in 'fixed' mode for a defined number of timesteps.
	sim_mode='fixed'
	
	if(length(metacomm)>0){
		comm = metacomm$comm
		poolA = metacomm$poolA
		poolB = metacomm$poolB
		sites = metacomm$sites
		assoc_probs = metacomm$assoc_probs
		gsad_a = metacomm$gsad_a
		gsad_b = metacomm$gsad_b
		niches_a = metacomm$niches_a
		niches_b = metacomm$niches_b
		topo = metacomm$topo
		topo_names = metacomm$topo_names
	}

	if(sim_mode=='fixed'){
	
	# Make vector of steps to save that includes timepoints in reps (if it is a vector)
	if(is.na(save_steps)) save_steps = 0:max(reps)
	save_steps = unique(c(save_steps, reps))
	save_steps = save_steps[order(save_steps)]

	# Generate empty array to hold changes in community and add initial community
	nsteps = length(save_steps)
	
	comm_records = array(NA, c(N_C, N, nsteps))
	poolA_records = array(NA, c(N_C, S_a, nsteps)) 
	poolB_records = array(NA, c(N_C, S_b, nsteps))
	Tmat_records = array(NA, c(N_C, N_L+1, N_L+1, nsteps))

	if(is.na(save_steps)|(0 %in% save_steps)){
		comm_records[,,1] = comm
		poolA_records[,,1] = poolA
		poolB_records[,,1] = poolB
	}

	# Run simulation
	for(step in 1:reps){
		# Mutualist mortality: mutualists present as associations cannot die
		poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
		poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)

		# Mutualists disperse into communities independently and with probability based on niche-based filtering	
		poolA = disperse(sites, niches_a, poolA, gsad_a)
		poolB = disperse(sites, niches_b, poolB, gsad_b)
	
		# Calculate transition matrices for each site
		T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs, omega) # array of square matrices of transition probabilities from one association to another

		# Identify current state of each space and transition based on random numbers
		new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
		for(i in 1:nrow(comm)){
		for(j in 1:ncol(comm)){
			new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
		}}
		comm = new_comm

		# Save communities
		if(step %in% save_steps){
			i = which(save_steps==step)
			comm_records[,,i] = comm
			poolA_records[,,i] = poolA
			poolB_records[,,i] = poolB
			Tmat_records[,,,i] = T_mat
		}
	}
	} # CLOSES if(sim_mode=='fixed') 

	# Run simulation until N chains converge on equilibrial community dynamics
	# Currently uses Rhat statistic from Brooks and Gelman 1998
	# This compares the within chain percentile interval to the interval calculated across all chains
#	if(sim_mode=='converge'){
#
#	# Split community into nchain identical instances
#	comm_arr = array(comm, dim=c(nrow(comm), ncol(comm), nchains))
#	poolA_arr = array(poolA, dim=c(nrow(poolA), ncol(poolA), nchains))
#	poolB_arr = array(poolB, dim=c(nrow(poolB), ncol(poolB), nchains))
#	Tmat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs)
#	Tmat_arr = array(Tmat, dim=c(dim(Tmat), nchains))
#	Rhat_arr = array()
#
#	# Generate empty array to hold changes in community and add initial community
#	comm_records = array(comm_arr, dim=c(dim(comm_arr), 1)) 
#	poolA_records = array(poolA_arr, dim=c(dim(poolA_arr),1)) 
#	poolB_records = array(poolB_arr, dim=c(dim(poolB_arr),1))
#	Tmat_records = array(Tmat_arr, dim=c(dim(Tmat_arr),1))
#
#	# Run simulation until convergence
#	window_size = 50 # Number of iterations to run before re-calculating convergence criteria
#	burnin = 100
#	thin = 1
#	Rhat_tol = 0.1 # Tolerance for R_hat statistic to be considered converged
#
#	# Set initial counters
#	converged = F
#	step = 1
#
#	# COMMUNITIES ARE CONVERGING VERY SLOWLY UNDER THESE CRITERIA- MAY WANT TO RECONSIDER
#	while(!converged){
#		# Do one step through the simulation
#		for(k in 1:nchains){
#			comm = comm_arr[,,k]
#			poolA = poolA_arr[,,k]
#			poolB = poolB_arr[,,k]
#
#			# Mutualist mortality: mutualists present as associations cannot die
#			poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
#			poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)
#
#			# Mutualists disperse into communities independently and with probability based on niche-based filtering	
#			poolA = disperse(sites, niches_a, poolA, gsad_a)
#			poolB = disperse(sites, niches_b, poolB, gsad_b)
#
#			# Calculate transition matrices for each site
#			T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs, omega) # array of square matrices of transition probabilities from one association to another
#			Tmat_arr[,,,k] = T_mat
#
#			# Identify current state of each space and transition based on random numbers
#			new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
#			for(i in 1:nrow(comm)){
#			for(j in 1:ncol(comm)){
#				new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
#			}} 
#			comm = new_comm
#
#			# Save communities
#			comm_arr[,,k] = comm
#			poolA_arr[,,k] = poolA
#			poolB_arr[,,k] = poolB
#		}
#
#		# Save communities from this step
#		Tmat_records = abind(Tmat_records, Tmat_arr)
#		comm_records = abind(comm_records, comm_arr)
#		poolA_records = abind(poolA_records, poolA_arr)
#		poolB_records = abind(poolB_records, poolB_arr)
#
#		# Calculate convergence criteria if past burnin and at correct re-test interval
#		if((step %% window_size == 0)&(step > burnin)){
#			# Determine observations to use in calculation of R_hat
#			window = floor(length((burnin+1):step)/2)
#			use_obs = seq((step-window+1), step, thin)
#			n_obs = length(use_obs)
#
#			# Calculate community composition~env correlation for each chain
#			corrstats = array(NA, dim=c(n_obs, nchains, 3, ncol(sites)), 
#				dimnames=list(step=use_obs, chain=1:nchains, community=c('species','a','b'), env=1:ncol(sites)))
#			for(n in use_obs){
#			for(k in 1:nchains){
#				corrstats[as.character(n), k, , ] = calc_envcorr(comm_records[,,k,n], topo_names, sites, c(), binary=F)
#			}}
#
#			# Calculate R_hat among chains
#			Rhats = sapply(dimnames(corrstats)$community, function(y){
#				sapply(1:ncol(sites), function(x){
#					calc_Rhat(corrstats[,,y,x], alpha=.95)
#				})
#			})
#
#			# Save records of Rhat
#			if(exists('Rhat_records')){ 
#				Rhat_records = abind(Rhat_records, Rhats)
#			} else {
#				Rhat_records = array(Rhats, dim=c(dim(Rhats),1), 
#					dimnames=list(env=1:dim(Rhats)[1], community=dimnames(Rhats)[[2]]))
#			}
#
#			corr_flags = abs(Rhats-1) < Rhat_tol
#
#			# Test whether all flags are TRUE
#			converged = prod(corr_flags)>0
#		}
#		step = step + 1

#	} # CLOSES while(converegence) loop
#	} # CLOSES if(sim_mode=='converge')

	# Return simulation results in a list
#	if(sim_mode=='converge') list(comm=comm_records, poolA=poolA_records, poolB=poolB_records, Tmat=Tmat_records, Rhat=Rhat_records)
	if(sim_mode=='fixed') list(comm=comm_records, poolA=poolA_records, poolB=poolB_records, Tmat=Tmat_records)
}
