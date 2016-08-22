

setwd('C:/Users/jrcoyle/Documents/Research/CAMM/')

results_dir = './Runs/Results_2/'
sim_dir = './GitHub/CAMM/R/'


# Load functions
source(paste0(sim_dir,'analysis_functions.R'))
source(paste0(sim_dir,'simulation_functions.R'))
source(paste0(sim_dir,'control_functions.R'))


# Examine a run
this_run = 'o-0_topo-many2many_envfilt-none_1'

load(paste0(results_dir, this_run, '_results.RData'))

topo_names = metacomm$topo_names
sites = metacomm$sites

# Calculate stats
rich_summary = sapply(1:nchains, function(i) calc_commstats(end_metacomms['comm',i][[1]], topo_names, list(a=end_metacomms['poolA',i][[1]], b=end_metacomms['poolB',i][[1]])), simplify='array')
comm_means = sapply(rownames(rich_summary), function(type) apply(simplify2array(rich_summary[type,]), 2, mean))
rich_stats = apply(comm_means, 2, function(x) c(mean(x), var(x))) 

rownames(rich_stats) = c('mean','var')



