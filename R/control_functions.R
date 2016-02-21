## This script holds functions to control the community assembly of mutualists model

# A function that writes a parameter file from a named list of parameters to a given directory
# Only works for parameters that are vectors of numerics or characters, lists with one level, and functions
# parm_list = a named list of parameters to write
# file_name = a string for the parameter filename
# file_dir = a string for the directory where to write the file
write_parms = function(parm_list, file_name, file_dir){
	make_newline = function(varname, value){
		if(length(value)==1){
			if(is.numeric(value)) new_line = paste(varname, value, sep='=')
			if(is.character(value)) new_line = paste0(varname, "='", value, "'")
		} else {
			if(is.numeric(value)) new_line = paste0(varname, '=c(', paste(value, collapse=','), ')')
			if(is.character(value)) new_line = paste0(varname, '=c(', paste0("'", value, "'", collapse=','), ')')
		}
		new_line
	}

	this_file = file(paste0(file_dir, file_name, '.txt'), open='w')
	for(i in 1:length(parm_list)){
		varname = names(parm_list)[i]
		value = parm_list[[i]]
	
		if(class(value)=='function'){
			new_line = c(paste0(varname,'='), deparse(value))	
		} else {	

		if(class(value)=='list'){
			line_list = paste(sapply(1:length(value), function(j){
				subvalue = value[[j]]
				make_newline(names(value)[j], subvalue)
			}), collapse=',')
			new_line = paste0(varname, '=list(', line_list, ')')
		} else {
			new_line = make_newline(varname, value)
		}}
				
		writeLines(new_line, this_file)
	}

	# Create lists of niche parameters
	writeLines('nicheparms_a = list(mu = c(mu_a1, mu_a2), rho = rho_a, sigma = c(sigma_a1, sigma_a2), alpha = c(alpha_a1, alpha_a2), r=r_a)', this_file)
	writeLines('nicheparms_b = list(mu = c(mu_b1, mu_b2), rho = rho_b, sigma = c(sigma_b1, sigma_b2), alpha = c(alpha_b1, alpha_b2), r=r_b)', this_file)

	close(this_file)	
}

# A function to make a list of parameters for writing to a file from a given environment
# This function needs to be updated whenever new parameters are coded into the simulation
make_parmlist = function(e=parent.frame()){
	list(
		runID = e$runID,
		a_name = e$a_name, 
		b_name = e$b_name,
		S_a = e$S_a,
		S_b = e$S_b,
		N_C = e$N_C,
		N = e$N,
		rho_z = e$rho_z, 
		mu_a1 = e$mu_a1,
		mu_a2 = e$mu_a2,
		rho_a = e$rho_a,
		mu_b1 = e$mu_b1,
		mu_b2 = e$mu_b2,
		rho_b = e$rho_b, 
		sigma_a1 = e$sigma_a1,
		sigma_a2 = e$sigma_a2,
		alpha_a1 = e$alpha_a1,
		alpha_a2 = e$alpha_a2, 
		r_a = e$r_a,
		sigma_b1 = e$sigma_b1, 
		sigma_b2 = e$sigma_b2,
		alpha_b1 = e$alpha_b1,
		alpha_b2 = e$alpha_b2,
		r_b = e$r_b,	
		gsad_dist_a = e$gsad_dist_a,
		gsad_dist_b = e$gsad_dist_b,
		gsad_cond_a = e$gsad_cond_a,
		gsad_cond_b = e$gsad_cond_b,
		N_L = e$N_L,
		topology = e$topology,
		omega = e$omega,
		assoc_prob_func = e$assoc_prob_func,
		mort_rate = e$mort_rate,
		mort_rate_a = e$mort_rate_a,
		mort_rate_b = e$mort_rate_b
	)
}




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
run_camm_N = function(sim_dir, parm_file, nruns, nchains, nparallel=1, sim_parms, simID, save_start=F, save_sim=F, save_dir='./'){
	runs = paste(simID, 1:nruns, sep='_')

	if(nparallel > 1){
		require(parallel)
		cluster = makeCluster(nparallel)	
		
		# Send required functions to each node
		clusterExport(cluster, c('runs','nchains','sim_dir','parm_file','sim_parms','simID','save_start','save_sim','save_dir'), envir=environment())
		clusterEvalQ(cluster, source(paste0(sim_dir, 'simulation_functions.R')))
		clusterEvalQ(cluster, source(parm_file))
		
		# Initialize CAMM
		metacomm_N = parLapply(cluster, 1:nruns, function(j) initialize_camm(parm_file, save_start, runID=runs[j], save_dir))
		clusterExport(cluster, 'metacomm_N', envir=environment())		
		
		# Run and Summarize CAMM
		sim_results = parLapply(cluster, 1:nruns, function(j){
			print(paste('start', runs[j]))
			metacomm = metacomm_N[[j]]			

			end_metacomms = sapply(1:nchains, function(i){
				# Run CAMM
				this_run = run_camm(metacomm=metacomm, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps[length(sim_parms$reps)], nchains=nchains)

				# Just output the instance of each community at each of the timepoints in reps
				list(comm=this_run$comm[,,sim_parms$reps+1], 
					poolA=this_run$poolA[,,sim_parms$reps+1],
					poolB=this_run$poolB[,,sim_parms$reps+1])
	
			})
			print(paste('done', runs[j]))

			# Define objects
			topo_names = metacomm$topo_names
			sites = metacomm$sites
	
			# Calculate mean community richness and abundance
			if(length(sim_parms$reps)==1){
				rich_summary = sapply(end_metacomms['comm',], function(comm) calc_commstats(comm, topo_names), simplify='array')
				comm_means = sapply(rownames(rich_summary), function(type) apply(simplify2array(rich_summary[type,]), 2, mean))
				rich_stats = apply(comm_means, 2, function(x) c(mean(x), var(x))) 
				rownames(rich_stats) = c('mean','var')
			
				# Calculate environmental correlations
				# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
				corr_stats = sapply(end_metacomms['comm',], function(comm) calc_envcorr(comm, topo_names, sites, 'jaccard', binary=T), simplify='array')
				corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
				dimnames(corr_stats)[[1]] = c('mean','var')
				# Returns array with [mean/var, type, env, S/N/rda/jaccard, binary]
			} else {
				rich_summary = sapply(end_metacomms['comm',], function(comms){
					 apply(comms, 3, function(comm) calc_commstats(comm, topo_names))	
				}, simplify='array')
				comm_means = apply(rich_summary, c(1,2), function(x) colMeans(x[[1]]))
				rich_stats = apply(comm_means, c(1,2), function(x) c(mean(x), var(x)))
				dimnames(rich_stats)[[1]] = c('mean','var')
				dimnames(rich_stats)[[3]] = paste0('T',sim_parms$reps)
				
				corr_stats = sapply(1:length(end_metacomms['comm',]), function(j){
					comms = end_metacomms['comm',][[j]]
					sapply(1:dim(comms)[3], function(i) calc_envcorr(comms[,,i], topo_names, sites, 'jaccard', binary=T), simplify='array')	
				}, simplify='array')
				corr_stats = apply(corr_stats, 1:5, function(vals) c(mean(vals), var(vals)))	
				dimnames(corr_stats)[[1]] = c('mean','var')
				dimnames(corr_stats)[[6]] = paste0('T',sim_parms$reps)
			}
			
			print(paste('summarized', runs[j]))
			# Save output
			if(save_sim) save(rich_stats, corr_stats, metacomm, end_metacomms, parm_file, sim_parms, nruns, nchains, file=paste0(save_dir, runs[j], '_results.RData'))
		
			list(rich_stats, corr_stats)	
		})
	
		stopCluster(cluster)

	} else {

		# Initialize CAMM
		metacomm_N = lapply(runs, function(x) initialize_camm(parm_file, save_start, runID=x, save_dir))
		
		# Run and Summarize CAMM
		sim_results = lapply(1:nruns, function(j){
			metacomm = metacomm_N[[j]]

			end_metacomms = sapply(1:nchains, function(i){
				# Run CAMM
				this_run = run_camm(metacomm, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps[length(sim_parms$reps)], nchains=sim_parms$nchains)
		
				# Just save the last instance of each community
				list(comm=this_run$comm[,,sim_parms$reps+1], 
					poolA=this_run$poolA[,,sim_parms$reps+1],
					poolB=this_run$poolB[,,sim_parms$reps+1])
			})

			# Define objects
			topo_names = metacomm$topo_names
			sites = metacomm$sites
	
			# Calculate mean community richness and abundance
			if(length(sim_parms$reps)==1){
				rich_summary = sapply(end_metacomms['comm',], function(comm) calc_commstats(comm, topo_names), simplify='array')
				comm_means = sapply(rownames(rich_summary), function(type) apply(simplify2array(rich_summary[type,]), 2, mean))
				rich_stats = apply(comm_means, 2, function(x) c(mean(x), var(x))) 
				rownames(rich_stats) = c('mean','var')
			
				# Calculate environmental correlations
				# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
				corr_stats = sapply(end_metacomms['comm',], function(comm) calc_envcorr(comm, topo_names, sites, 'jaccard', binary=T), simplify='array')
				corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
				dimnames(corr_stats)[[1]] = c('mean','var')
				# Returns array with [mean/var, type, env, S/N/rda/jaccard, binary]
			} else {
				rich_summary = sapply(end_metacomms['comm',], function(comms){
					 apply(comms, 3, function(comm) calc_commstats(comm, topo_names))	
				}, simplify='array')
				comm_means = apply(rich_summary, c(1,2), function(x) colMeans(x[[1]]))
				rich_stats = apply(comm_means, c(1,2), function(x) c(mean(x), var(x)))
				dimnames(rich_stats)[[1]] = c('mean','var')
				dimnames(rich_stats)[[3]] = paste0('T',sim_parms$reps)
				
				corr_stats = sapply(1:length(end_metacomms['comm',]), function(j){
					comms = end_metacomms['comm',][[j]]
					sapply(1:dim(comms)[3], function(i) calc_envcorr(comms[,,i], topo_names, sites, 'jaccard', binary=T), simplify='array')	
				}, simplify='array')
				corr_stats = apply(corr_stats, 1:5, function(vals) c(mean(vals), var(vals)))	
				dimnames(corr_stats)[[1]] = c('mean','var')
				dimnames(corr_stats)[[6]] = paste0('T',sim_parms$reps)
			}
			
			# Save output
			if(save_sim) save.image(paste0(save_dir, runs[j], '_results.RData'))

			list(rich_stats, corr_stats)	
		})
	}
	
	sim_results

} 


# A function that summarizes a given parameter for a simulation run with run_camm_N()
# Returns an array with the mean, variance, median, and 95th percentile intervals as the 1st dimension
# results = a list of simulation results returns by run_camm_N
# what = a string indicating which community chacteristic should be summarized
# type = a string indicating whether statistics should be computed for 'species', 'a', or 'b'
summarize_camm = function(results, what, type=NA){
	# Determine whether results have multiple time points
	tp = length(dim(results[[1]][[1]]))==3

	if(tp){
		# Richness in each community
		if(what=='S'){
			get_vars = paste(c('S','Stot'),type,sep='_')
			richness = sapply(results, function(x) x[[1]][,get_vars,], simplify='array')
			return_stats = apply(richness, 1:3, function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response','time')
		}
		
		if(what=='N'){
			abun = sapply(results, function(x) x[[1]][,c('N', 'Ntot'),], simplify='array')
			return_stats = apply(abun, 1:3, function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response','time')
		}
		
		if(what=='cor'){
			corr_mat = sapply(results, function(x) x[[2]][,type,,,,], simplify='array')
			maxD = length(dim(corr_mat))
			return_stats = apply(corr_mat, 1:(maxD-1), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','env','measure','time')
		}
	} else {
		# Richness in each community
		if(what=='S'){
			get_vars = paste(c('S','Stot'),type,sep='_')
			richness = sapply(results, function(x) x[[1]][,get_vars], simplify='array')
			return_stats = apply(richness, c(1,2), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response')
		}


		if(what=='N'){
			abun = sapply(results, function(x) x[[1]][,c('N', 'Ntot')], simplify='array')
			return_stats = apply(abun, c(1,2), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response')
		}

		if(what=='cor'){
			corr_mat = sapply(results, function(x) x[[2]][,type,,,], simplify='array')
			maxD = length(dim(corr_mat))
			return_stats = apply(corr_mat, 1:(maxD-1), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','env','measure')
		}
	}

	return_stats
}





