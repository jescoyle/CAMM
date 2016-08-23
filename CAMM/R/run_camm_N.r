#' Run multiple simulations
#'
#' Runs multiple simulations on the same set of parameters.
#'
#' This function runs multiple simulations on a set of parameters 
#' given in \code{parm_file}. The function first generates \code{nruns}
#' independent metacommunities and if \code{save_start} is \code{TRUE}, 
#' saves these as a list in the directory given in \code{save_dir} as
#' \code{'<simID>_metacomms.RData'}. The function then executes one or more
#' simulations (given in \code{nchains}) on each of these metacommunities.
#' This allows the user to evaluate variability among simulations with 
#' different initial metacommunities versus variability that arises during 
#' the simulation. Simulations are controled by parameters passed in 
#' \code{sim_parms}, including the number of timesteps and which steps to 
#' record. Mean diversity statistics across communities ('rich_stats') and  
#' community-environment correlations ('corr_stats') are calculated for  
#' each simulation using \code{\link{calc_commstats}} and 
#' \code{\link{calc_envcorr}}) and the mean and variance of these community 
#' descriptors are calculated across simulations run on the same initial  
#' metacommunity (across 'chains'). Summaries can be calculated for multiple  
#' timepoints by passing a vector of timepoints to the \code{sim_parms$reps}.   
#' In this case, the simulation will be run until the last timepoint in the  
#' vector (not the largest). If \code{save_sim} is \code{TRUE} then both the  
#' cross-chain summary statistics and full simulation results are saved to 
#' files in \code{save_dir} as \code{'<simID>_<run#>_results.RData'} and 
#' \code{'<simID>_<run#>_metacomms-end.RData'}, respectively. Otherwise
#' only the summaries are returned by the function, in a list.
#'
#' Simulations can be run in parallel by specifying 
#' \code{nparallel > 1}, which requires the \code{\link[doParallel]{doParallel}}
#' and \code{\link[foreach]{foreach}} packages.
#' By default, \code{nparallel = 1} and the simulations proceed serially.
#'
#' @param parm_file (required) path to file where parameters are stored. 
#' 	See \code{\link{write_parms}} and \code{\link{make_parmlist}} for 
#' 	instructions on creating an appropriate parameter file
#' @param nruns (required) number of independent simulation runs to conduct
#' @param nchains (required) number of replicate simulations on each initial
#' 	metacommunity. Defaults to 1.
#' @param nparallel number of cores used for parallel simulation. Defaults to
#' 1.
#' @param sim_parms (required) named list of parameters controlling simulation.
#' 	Should contain:
#' 	\describe{
#' 		\item{reps}{vector of timepoints at which to summarize simulations
#' 			across chains (see details). The last number is the timepoint at
#' 			which simulations are stopped.}
#' 		\item{save_steps}{optional vector of timepoints to save. Defaults to 
#' 			all.}
#' 	}
#' @param simID (required) character vector identifying this simulation
#' @param save_start logical indicating whether the initial metacommunities
#' 	should be saved to an RData file. Defaults to \code{FALSE}.
#' @param save_sim logical indicating whether simulations results should
#' 	be save to RData files. Defaults to \code{FALSE}.
#' @param save_dir path to directory where simulation results should be saved.
#' 	Defaults to the current directory.
#' @param restart logical indicating whether this simulation restarts an existing
#' 	simulation. If \code{FALSE}, existing files will be overwritten. If 
#' 	\code{TRUE}, then a metacommunity object named 
#' 	\code{'<runID>_metacomms.RData'} must be present in \code{save_dir}.
#' 	Defaults to \code{FALSE}.
#' @param sim_dir directory where CAMM is installed. Defaults to R's search path.
#' @param dev logical flag indicating whether the function is being called in
#' 	development mode, in which case the package is loaded using devtools from
#' 	\code{sim_dir}
#' @return a list of lists (one for each independent simulation run), each 
#' 	containing two arrays summarizing diversity ('rich_stats') and 
#' 	community-environment correlations ('corr_stats') across simulations run   
#' 	on the same initial metacommunity (across chains). The first dimension of 
#' 	these arrrays refers to whether the quantities are the cross-chain mean
#' 	or variance. If summaries were requested for multiple timepoints, then 
#' 	the timepoint is indicated in the last dimension. 
#'
#' @export

run_camm_N = function(parm_file, nruns, nchains, nparallel=1, sim_parms, simID, save_start=F, save_sim=F, save_dir='./', restart=F, sim_dir=NULL, dev=FALSE){
	runs = paste(simID, 1:nruns, sep='_')

	# Make directory for saving results
	if(!file.exists(save_dir)) dir.create(save_dir)
	
	# Attempt to run simulations in parallel
	if((nparallel > 1) & requireNamespace('parallel', quietly=TRUE)){
		# Attach functions in parallel
		#library(parallel)

		cluster = parallel::makeCluster(nparallel, outfile=file.path(save_dir, paste0(simID, '.Rout')))	
		
		# Send required objects and functions to each node
		parallel::clusterExport(cluster, c('runs','nchains','sim_dir','parm_file','sim_parms','simID','save_start','save_sim','save_dir','restart'), envir=environment())

		# During development, must load development version of package
		# Otherwise, load package from sim_dir or default R search path
		if(dev){
			parallel::clusterEvalQ(cluster, {
				library(devtools)
				current_code = as.package(file.path(sim_dir,'CAMM'))
				load_all(current_code)
			})
		} else {
			parallel::clusterEvalQ(cluster, library(CAMM, lib.loc=sim_dir))
		}

		parallel::clusterEvalQ(cluster, source(parm_file))

		# Initialize CAMM or get previously saved initial metacommunities
		if(restart){
			# Get pre-existing metacomm_N object
			load(file.path(save_dir, paste0(simID, '_metacomms.RData')))
		} else {
			#print('starting runs')
			metacomm_N = parallel::parLapply(cluster, 1:nruns, function(j) initialize_camm(parm_file, save_start, runID=runs[j], save_dir))
		}
		parallel::clusterExport(cluster, 'metacomm_N', envir=environment())

		if(save_start&(!restart)) save(metacomm_N, file=file.path(save_dir, paste0(simID, '_metacomms.RData')))

		# Run and Summarize CAMM
		sim_results = parallel::parLapply(cluster, 1:nruns, function(j){
		
			# Load parameters into this environment
			#source(parm_file)
			
			print(paste('start', runs[j]))
			print(environment())
			print(parent.frame())
			metacomm = metacomm_N[[j]]
			reps = sim_parms$reps

			# Define steps to save during simulation that includes timepoints in reps
			if(is.null(sim_parms$save_steps)){
				save_steps = 0:max(reps)
			} else {
				save_steps = sim_parms$save_steps
			}
			save_steps = unique(c(save_steps, reps))
			save_steps = save_steps[order(save_steps)]
	
			# Check whether run exists and load results if it does and this is a restart
			# Otherwise, do run.
			this_file = file.path(save_dir, paste0(runs[j], '_metacomms-end.RData'))
			if(restart & file.exists(this_file)){
				load(this_file)
			} else {
				
				end_metacomms = sapply(1:nchains, function(i){
					print(paste('starting chain',i))
					print(environment())
					print(parent.frame())
					print(parent.frame(2))
					print(parent.frame(3))
					
					# Run CAMM
					this_run = run_camm(metacomm=metacomm, reps=sim_parms$reps[length(sim_parms$reps)], save_steps=save_steps)

					# Just output the instance of each community at each of the timepoints in reps
					reps_i = which(save_steps %in% reps)
					list(comm=this_run$comm[,,reps_i], 
						poolA=this_run$poolA[,,reps_i],
						poolB=this_run$poolB[,,reps_i])
				})
				if(save_sim) save(end_metacomms, file=this_file)
			}

			# Check whether summary exists load summary if it does and this is a restart
			this_file = file.path(save_dir, paste0(runs[j], '_results.RData'))
			if(restart & file.exists(this_file)){
				load(this_file)
			} else {
				# Define objects
				topo_names = metacomm$topo_names
				sites = metacomm$sites
		
				# Calculate mean community richness and abundance
				# For runs where only the end timepoint should be summarized
				if(length(sim_parms$reps)==1){
					rich_summary = sapply(1:nchains, function(i) calc_commstats(end_metacomms['comm',i][[1]], topo_names, list(a=end_metacomms['poolA',i][[1]], b=end_metacomms['poolB',i][[1]])), simplify='array')
					comm_means = sapply(rownames(rich_summary), function(type) apply(simplify2array(rich_summary[type,]), 2, mean))
					rich_stats = apply(comm_means, 2, function(x) c(mean(x), var(x))) 
					rownames(rich_stats) = c('mean','var')
				
					# Calculate environmental correlations
					# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
					corr_stats = sapply(1:nchains, function(i) calc_envcorr(end_metacomms['comm',i][[1]], topo_names, sites, 'jaccard', binary=T,list(a=end_metacomms['poolA',i][[1]], b=end_metacomms['poolB',i][[1]])), simplify='array')
					corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
					dimnames(corr_stats)[[1]] = c('mean','var')
					# Returns array with [mean/var, type, env, S/N/rda/jaccard, binary]

				# For runs where multiple timepoints should be summarized
				} else {
					rich_summary = sapply(1:nchains, function(i){
						comms = end_metacomms[,i]$comm
						poolAs = end_metacomms[,i]$poolA
						poolBs = end_metacomms[,i]$poolB
						sapply(1:dim(comms)[3], function(j) calc_commstats(comms[,,j], topo_names, list(a=poolAs[,,j], b=poolBs[,,j])))	
					}, simplify='array')
					comm_means = apply(rich_summary, c(1,2,3), function(x) mean(x[[1]]))
					rich_stats = apply(comm_means, c(1,2), function(x) c(mean(x), var(x)))
					dimnames(rich_stats)[[1]] = c('mean','var')
					dimnames(rich_stats)[[3]] = paste0('T',sim_parms$reps)
					
					corr_stats = sapply(1:nchains, function(i){
						comms = end_metacomms[,i]$comm
						poolAs = end_metacomms[,i]$poolA
						poolBs = end_metacomms[,i]$poolB
						sapply(1:dim(comms)[3], function(j) calc_envcorr(comms[,,j], topo_names, sites, 'jaccard', binary=T, list(a=poolAs[,,j], b=poolBs[,,j])), simplify='array')	
					}, simplify='array')
					corr_stats = apply(corr_stats, 1:5, function(vals) c(mean(vals), var(vals)))	
					dimnames(corr_stats)[[1]] = c('mean','var')
					dimnames(corr_stats)[[6]] = paste0('T',sim_parms$reps)
				}

				print(paste('summarized', runs[j]))
				# Save output
				if(save_sim) save(rich_stats, corr_stats, metacomm, end_metacomms, parm_file, sim_parms, nruns, nchains, file=this_file)
			}
			
			list(rich_stats, corr_stats)
		})

		parallel::stopCluster(cluster)

	} else {
	
		if(nparallel>1) warning('parallel package not found. Running simulations serially.')

		# Initialize CAMM or get previously saved initial metacommunities
		if(restart){
			# Get pre-existing metacomm_N object
			load(file.path(save_dir, paste0(simID, '_metacomms.RData')))
		} else {
			metacomm_N = lapply(runs, function(x) initialize_camm(parm_file, save_start, runID=x, save_dir))
		}
		
		if(save_start&(!restart)) save(metacomm_N, file=file.path(save_dir, paste0(simID, '_metacomms.RData')))
		
		# Run and Summarize CAMM
		sim_results = lapply(1:nruns, function(j){
			metacomm = metacomm_N[[j]]
			reps = sim_parms$reps

			# Define steps to save during simulation that includes timepoints in reps
			if(is.null(sim_parms$save_steps)){
				save_steps = 0:max(reps)
			} else {
				save_steps = sim_parms$save_steps
			}
			save_steps = unique(c(save_steps, reps))
			save_steps = save_steps[order(save_steps)]
	
			# Check whether run exists and load results if it does and this is a restart
			# Otherwise, do run.
			this_file = file.path(save_dir, paste0(runs[j], '_metacomms-end.RData'))
			if(restart & file.exists(this_file)){
				load(this_file)
			} else {
				end_metacomms = sapply(1:nchains, function(i){
					# Run CAMM
					this_run = run_camm(metacomm=metacomm, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps[length(sim_parms$reps)], 
						nchains=nchains, save_steps=save_steps)

					# Just output the instance of each community at each of the timepoints in reps
					reps_i = which(save_steps %in% reps)
					list(comm=this_run$comm[,,reps_i], 
						poolA=this_run$poolA[,,reps_i],
						poolB=this_run$poolB[,,reps_i])
		
				})
				if(save_sim) save(end_metacomms, file=this_file)
			}

			# Check whether summary exists load summary if it does and this is a restart
			this_file = file.path(save_dir, paste0(runs[j], '_results.RData'))
			if(restart & file.exists(this_file)){
				load(this_file)
			} else {

				# Define objects
				topo_names = metacomm$topo_names
				sites = metacomm$sites
		
				# Calculate mean community richness and abundance
				if(length(sim_parms$reps)==1){
					rich_summary = sapply(1:nchains, function(i) calc_commstats(end_metacomms['comm',i][[1]], topo_names, list(a=end_metacomms['poolA',i][[1]], b=end_metacomms['poolB',i][[1]])), simplify='array')
					comm_means = sapply(rownames(rich_summary), function(type) apply(simplify2array(rich_summary[type,]), 2, mean))
					rich_stats = apply(comm_means, 2, function(x) c(mean(x), var(x))) 
					rownames(rich_stats) = c('mean','var')
				
					# Calculate environmental correlations
					# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
					corr_stats = sapply(1:nchains, function(i) calc_envcorr(end_metacomms['comm',i][[1]], topo_names, sites, 'jaccard', binary=T,list(a=end_metacomms['poolA',i][[1]], b=end_metacomms['poolB',i][[1]])), simplify='array')
					corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
					dimnames(corr_stats)[[1]] = c('mean','var')
					# Returns array with [mean/var, type, env, S/N/rda/jaccard, binary]
				} else {
					rich_summary = sapply(1:nchains, function(i){
						comms = end_metacomms[,i]$comm
						poolAs = end_metacomms[,i]$poolA
						poolBs = end_metacomms[,i]$poolB
						sapply(1:dim(comms)[3], function(j) calc_commstats(comms[,,j], topo_names, list(a=poolAs[,,j], b=poolBs[,,j])))	
					}, simplify='array')
					comm_means = apply(rich_summary, c(1,2,3), function(x) mean(x[[1]]))
					rich_stats = apply(comm_means, c(1,2), function(x) c(mean(x), var(x)))
					dimnames(rich_stats)[[1]] = c('mean','var')
					dimnames(rich_stats)[[3]] = paste0('T',sim_parms$reps)
					
					corr_stats = sapply(1:nchains, function(i){
						comms = end_metacomms[,i]$comm
						poolAs = end_metacomms[,i]$poolA
						poolBs = end_metacomms[,i]$poolB
						sapply(1:dim(comms)[3], function(j) calc_envcorr(comms[,,j], topo_names, sites, 'jaccard', binary=T, list(a=poolAs[,,j], b=poolBs[,,j])), simplify='array')	
					}, simplify='array')
					corr_stats = apply(corr_stats, 1:5, function(vals) c(mean(vals), var(vals)))	
					dimnames(corr_stats)[[1]] = c('mean','var')
					dimnames(corr_stats)[[6]] = paste0('T',sim_parms$reps)
				}
				
				# Save output
				if(save_sim) save(rich_stats, corr_stats, metacomm, end_metacomms, parm_file, sim_parms, nruns, nchains, file=this_file)
			}

			list(rich_stats, corr_stats)
		})
	}
	
	# Return simulations results
	sim_results
} 

