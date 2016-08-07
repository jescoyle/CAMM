### This script was used to re-make run summaries from early runs

options(stringsAsFactors=F)
library(reshape) # melt()
library(foreach)
library(doParallel)

# The following arguments should be passed (in order)
# directory where old results  saved: defaults to './Results'
# directory with simulation scripts : defaults to current directory
# directory where new results should be saved : defaults to './Results'
args = commandArgs(trailingOnly=T)
print(args)

# Set directories
results_dir = ifelse(is.na(args[1]), './Results/', args[1])
sim_dir = ifelse(is.na(args[2]), './', args[2])
save_dir = ifelse(is.na(args[3]), './Results/', args[3])

# Set number of cores
ncores = ifelse(is.na(args[4]), 1, as.numeric(args[4]))

# Load functions
source(paste0(sim_dir,'analysis_functions.R'))
source(paste0(sim_dir,'simulation_functions.R'))
source(paste0(sim_dir,'control_functions.R'))

# Get list of results files
results_filelist = list.files(results_dir, 'results.RData')

# Find unique parameter combinations
runIDs = sapply(results_filelist, function(x) gsub('_[0-9]+_results.RData', '', x))
runIDs = unique(runIDs)

# Make cluster
cluster = makeCluster(ncores)
registerDoParallel(cluster)

foreach(runID=runIDs, .packages=c('reshape','vegan')) %dopar% {
	# Get list of results files from this set of parameters
	filelist = list.files(results_dir, paste0(runID,'_[0-9]+_results.RData'))

	# Extract parameter values
	parm_vals = get_parms(runID)

	model_out = lapply(filelist, function(f){
		load(paste0(results_dir, f))

		# Define objects that define this particular simulation
		topo_names = metacomm$topo_names
		sites = metacomm$sites

		# Calculate mean community richness and abundance
		# For runs where only the end tiepoint shoudl be summarized
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
		
		list(rich_stats, corr_stats)
	})

	# Analyze results
	S_a = melt(summarize_camm(model_out, 'S', 'a'))
	S_b = melt(summarize_camm(model_out, 'S', 'b'))
	N_comm = melt(summarize_camm(model_out, 'N'))
	Cor_ab = melt(summarize_camm(model_out, 'Cor_ab'))
	Cor_ab$response = 'Cor_ab'
	Cor_ab = Cor_ab[,colnames(S_a)]
	comm_summary = rbind(S_a, S_b, Cor_ab, N_comm)
	comm_summary = cast(comm_summary, ... ~ response)
	
	cor_a = melt(summarize_camm(model_out, 'cor','a'))
	names(cor_a)[names(cor_a)=='value'] = 'cor_a'
	cor_b = melt(summarize_camm(model_out, 'cor','b'))
	names(cor_b)[names(cor_b)=='value'] = 'cor_b'
	cor_summary = merge(cor_a, cor_b)

	write.csv(comm_summary, file=paste0(save_dir, runID, '_comm_summary.csv'), row.names=F)
	write.csv(cor_summary, file=paste0(save_dir, runID, '_cor_summary.csv'), row.names=F)

	i_run = which(runIDs==runID)
	print(paste('Finished', i_run, '/', length(runIDs)))
}

stopCluster(cluster)


