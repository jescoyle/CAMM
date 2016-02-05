## This script holds functions to control the community assembly of mutualists model


# A function that starts multiple runs of a simulation
# sim_dir = a directory containing simulations_functions.R
# parm_file = a file containing parameter values
# nruns = number of simulation runs to start
# nchains = number of identical starts in each run
# nparallel = if run in parallel, number of cores
# sim_parms = named list of parameter controling when simulations should stop
#	components are: sim_mode and reps or nchains (depending on sim_mode)
# simID = string that identifies this set of simulations runs
# save_start = flag indicating whether each initial set of niches and topologies should be saved
run_camm_N = function(sim_dir, parm_file, nruns, nchains, nparallel=1, sim_parms, simID, save_start=F){
	runs = paste0(simID, 1:nruns)

	if(nparallel > 1){
		require(parallel)
		cluster = makeCluster(nparallel)	
	
		# Send required functions to each node
		clusterExport(cluster, c('runs','nchains','sim_dir','parm_file','sim_parms','simID','save_start'))
		clusterEvalQ(cluster, source(paste0(sim_dir, 'simulation_functions.R')))
		clusterEvalQ(cluster, source(parm_file))

		# Initialize CAMM
		metacomm_N = parLapply(cluster, runs, function(x) initialize_camm(parm_file, save_start, runID=x))
		
		# Run and Summarize CAMM
		sim_results = parLapply(cluster, metacomm_N, function(x){
			end_metacomms = sapply(1:nchains, function(i){
				# Run CAMM
				this_run = run_camm(metacomm=x, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps, nchains=sim_parms$nchains)
		
				# Just save the last instance of each community
				list(comm=this_run$comm[,,sim_parms$reps+1], 
					poolA=this_run$poolA[,,sim_parms$reps+1],
					poolB=this_run$poolB[,,sim_parms$reps+1])
			})
			
			# Define objects
			topo_names = x$topo_names
			sites = x$sites
	
			# Calculate richness and abundance
			rich_stats = sapply(end_metacomms['comm',], function(comm) calc_commstats(comm, topo_names), simplify='array')
			rich_stats = sapply(rownames(rich_stats), function(type) apply(simplify2array(rich_stats[type,]), 1, function(x) c(mean(x), var(x))), simplify='array') 
			dimnames(rich_stats)[[1]] = c('mean','var')
			dimnames(rich_stats)[[2]] = 1:N_C
			# Returns array with [mean/var, community, richness/abundance]

			# Calculate environmental correlations
			# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
			corr_stats = sapply(end_metacomms['comm',], function(comm) calc_envcorr(comm, topo_names, sites, 'jaccard', binary=T), simplify='array')
			corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
			dimnames(corr_stats)[[1]] = c('mean','var')
			# Returns array with [mean/var, type, env, rda/jaccard, binary]

			list(rich_stats, corr_stats)	
		})
	
		stopCluster(cluster)

	} else {

		# Initialize CAMM
		metacomm_N = lapply(runs, function(x) initialize_camm(parm_file, save_start, runID=x))
		
		# Run and Summarize CAMM
		sim_results = lapply(metacomm_N, function(x){
			end_metacomms = sapply(1:nchains, function(i){
				# Run CAMM
				this_run = run_camm(metacomm=x, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps, nchains=sim_parms$nchains)
		
				# Just save the last instance of each community
				list(comm=this_run$comm[,,sim_parms$reps+1], 
					poolA=this_run$poolA[,,sim_parms$reps+1],
					poolB=this_run$poolB[,,sim_parms$reps+1])
			})
			
			# Define objects
			topo_names = x$topo_names
			sites = x$sites
	
			# Calculate richness and abundance
			rich_stats = sapply(end_metacomms['comm',], function(comm) calc_commstats(comm, topo_names), simplify='array')
			rich_stats = sapply(rownames(rich_stats), function(type) apply(simplify2array(rich_stats[type,]), 1, function(x) c(mean(x), var(x))), simplify='array') 
			dimnames(rich_stats)[[1]] = c('mean','var')
			dimnames(rich_stats)[[2]] = 1:N_C
			# Returns array with [mean/var, community, richness/abundance]

			# Calculate environmental correlations
			# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
			corr_stats = sapply(end_metacomms['comm',], function(comm) calc_envcorr(comm, topo_names, sites, 'jaccard', binary=T), simplify='array')
			corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
			dimnames(corr_stats)[[1]] = c('mean','var')
			# Returns array with [mean/var, type, env, rda/jaccard, binary]

			list(rich_stats, corr_stats)	
		})
	}

} 



