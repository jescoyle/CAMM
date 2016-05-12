### This script was used to re-make run summaries from early runs

options(stringsAsFactors=F)
library(reshape) # melt()

# The following arguments should be passed (in order)
# directory where results  saved: defaults to './Results'
# directory with simulation scripts : defaults to current directory
# directory where new results should be saved : defaults to './Results'
args = commandArgs(trailingOnly=T)
print(args)

# Set directories
results_dir = ifelse(is.na(args[1]), './Results/', args[1])
sim_dir = ifelse(is.na(args[2]), './', args[2])
save_dir = ifelse(is.na(args[3]), './Results/', args[3])


# Load functions
code_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/GitHub/CAMM/R/'
source(paste0(sim_dir,'analysis_functions.R'))
source(paste0(sim_dir,'simulation_functions.R'))
source(paste0(sim_dir,'control_functions.R'))

# Get list of results files
results_filelist = list.files(results_dir, 'results.RData')

# Find unique parameter combinations
runIDs = sapply(results_filelist, function(x) gsub('_[0-9]+_results.RData', '', x))
runIDs = unique(runIDs)

for(runID in runIDs){
	# Get list of results files from this set of parameters
	filelist = list.files(results_dir, paste0(runID,'_[0-9]+_results.RData'))

	# Extract parameter values
	parm_vals = get_parms(runID)

	model_out = lapply(filelist, function(f){
		load(paste0(results_dir, f))

		# Define objects that define this particular simulation
		topo_names = metacomm$topo_names
		sites = metacomm$sites

		# Summarize results across chains
		# This works if reps was a single integer (only one endpoint saved)
		rich_summary = sapply(end_metacomms['comm',], function(comm) calc_commstats(comm, topo_names), simplify='array')
		comm_means = sapply(rownames(rich_summary), function(type) apply(simplify2array(rich_summary[type,]), 2, mean))
		rich_stats = apply(comm_means, 2, function(x) c(mean(x), var(x))) 
		rownames(rich_stats) = c('mean','var')
			
		# Calculate environmental correlations
		# CURRENTLY USES PRES/ABSENCE, NOT ABUNDANCE
		corr_stats = sapply(end_metacomms['comm',], function(comm) calc_envcorr(comm, topo_names, sites, 'jaccard', binary=T), simplify='array')
		corr_stats = apply(corr_stats, 1:4, function(vals) c(mean(vals), var(vals)))	
		dimnames(corr_stats)[[1]] = c('mean','var')
		
		list(rich_stats, corr_stats)
	})

	model_out = sim_results

	# Analyze results
	S_a = melt(summarize_camm(model_out, 'S', 'a'))
	S_b = melt(summarize_camm(model_out, 'S', 'b'))
	N_comm = melt(summarize_camm(model_out, 'N'))
	comm_summary = rbind(S_a, S_b, N_comm)
	comm_summary = cast(comm_summary, ... ~ response)
	
	cor_a = melt(summarize_camm(model_out, 'cor','a'))
	names(cor_a)[names(cor_a)=='value'] = 'cor_a'
	cor_b = melt(summarize_camm(model_out, 'cor','b'))
	names(cor_b)[names(cor_b)=='value'] = 'cor_b'
	cor_summary = merge(cor_a, cor_b)

	write.csv(comm_summary, file=paste0(save_dir, runID, '_comm_summary.csv'), row.names=F)
	write.csv(cor_summary, file=paste0(save_dir, runID, '_cor_summary.csv'), row.names=F)
}


