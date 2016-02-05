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
run_camm_N = function(sim_dir, parm_file, nruns, nchains, nparallel=1, sim_parms, simID, save_start=F, save_start=F){
	runs = paste0(simID, 1:nruns)

	if(nparallel > 1){
		require(parallel)
		cluster = makeCluster(nparallel)	
	
		# Send required functions to each node
		clusterExport(cluster, c('runs','nchains','sim_dir','parm_file','sim_parms','simID','save_start','save_sim'))
		clusterEvalQ(cluster, source(paste0(sim_dir, 'simulation_functions.R')))
		clusterEvalQ(cluster, source(parm_file))

		# Initialize CAMM
		metacomm_N = parLapply(cluster, 1:nruns, function(j) initialize_camm(parm_file, save_start, runID=runs[j]))
		clusterExport(cluster, 'metacomm_N')		

		# Run and Summarize CAMM
		sim_results = parLapply(cluster, 1:nruns, function(j){
			metacomm = metacomm_N[[j]]			

			end_metacomms = sapply(1:nchains, function(i){
				# Run CAMM
				this_run = run_camm(metacomm, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps, nchains=sim_parms$nchains)

				# Just output the last instance of each community
				list(comm=this_run$comm[,,sim_parms$reps+1], 
					poolA=this_run$poolA[,,sim_parms$reps+1],
					poolB=this_run$poolB[,,sim_parms$reps+1])
			})

			# Save output
			if(save_sim) save.image(paste(runs[j], 'results.RData', sep='_'))
			
			# Define objects
			topo_names = metacomm$topo_names
			sites = metacomm$sites
	
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
		sim_results = lapply(1:nruns, function(j){
			metacomm = metacomm_N[[j]]

			end_metacomms = sapply(1:nchains, function(i){
				# Run CAMM
				this_run = run_camm(metacomm, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps, nchains=sim_parms$nchains)
		
				# Just save the last instance of each community
				list(comm=this_run$comm[,,sim_parms$reps+1], 
					poolA=this_run$poolA[,,sim_parms$reps+1],
					poolB=this_run$poolB[,,sim_parms$reps+1])
			})

			# Save output
			if(save_sim) save.image(paste(runs[j], 'results.RData', sep='_'))
			
			# Define objects
			topo_names = metacomm$topo_names
			sites = metacomm$sites
	
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
	
	sim_results

} 


# A function that summarizes a given parameter for a simulation run with run_camm_N()
# Returns an array with the mean, variance, median, and 95th percentile intervals as the 1st dimension
# results = a list of simulation results returns by run_camm_N
# what = a string indicating which community chacteristic should be summarized
# type = a strong indicating whether statistics should be computed for 'species', 'a', or 'b'
summarize_camm = function(results, what, type){
	# Richness in each community
	if(what=='S'){
		get_var = paste('S',type,sep='_')
		richness = sapply(results, function(x) x[[1]][,,get_var], simplify='array')
		return_stats = apply(richness, c(1,2), function(x) c(mean=mean(x), var=var(x), quantile(x, c(0.025, 0.5, 0.975))))
	}

	if(what=='N'){
		abun = sapply(results, function(x) x[[1]][,,'N'], simplify='array')
		return_stats = apply(abun, c(1,2), function(x) c(mean=mean(x), var=var(x), quantile(x, c(0.025, 0.5, 0.975))))
	}

	if(what=='cor'){
		corr_mat = sapply(results, function(x) x[[2]][,type,,,], simplify='array')
		maxD = length(dim(corr_mat))
		return_stats = apply(corr_mat, 1:(maxD-1), function(x) c(mean=mean(x), var=var(x), quantile(x, c(0.025, 0.5, 0.975))))
	}

	return_stats
}





