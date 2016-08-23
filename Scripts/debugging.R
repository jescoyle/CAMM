# This script is used for debugging CAMM functions

#setwd('C:/Users/jrcoyle/Documents/Research/CAMM/GitHub')
#install.packages('CAMM_0.1.0.zip', repos=NULL)

library(CAMM)

options(stringsAsFactors=F)
setwd('C:/Users/jrcoyle/Documents/Research/CAMM')

# Load development version
library(devtools)
library(roxygen2)
current_code = as.package('GitHub/CAMM')
load_all(current_code)


# Set directory for saving files
save_dir = 'Sandbox'

# Set simulation parameter file
parmfile = 'Sandbox/parameter_file.R'
parmfile = 'Sandbox/MP/NL-90_envfilt-none/parameter_file.R'

# Check initialize camm
mc = initialize_camm(parmfile, save_start = F, runID='test', save_dir='Sandbox')


# Set simulation controls
nruns = 4
nchains = 3
sim_parms = list(reps = c(5,10))
simID = 'test'


# Run simulation
results = run_camm_N(parmfile, nruns, nchains, nparallel=2, sim_parms, simID, 
	save_start=T, save_sim=T, save_dir, sim_dir='GitHub', dev=T)

# Example summary
summary = summarize_camm(results, 'S', 'a')


# Load metacommunity objects
load(file.path(save_dir, paste0(simID, '_metacomms.RData')))
load(file.path(save_dir, paste0(simID, '_1_metacomms-end.RData')))

metacomm_N[[1]]

# Assign objects to work with
comm = end_metacomms['comm',1]$comm[,,2] # 1st chain, last timepoint
topo = metacomm_N[[1]]$topo
topo_names = metacomm_N[[1]]$topo_names

calc_network(comm, topo_names, bysite=F)

netA = calc_network(comm, topo_names, bysite=T)[,,1]
netB = calc_network(comm, topo_names, bysite=T)[,,2]

compare_networks(netA, netB)

site_nets = calc_network(comm, topo_names, bysite=T)
full_net = calc_network(comm, topo_names)

sapply(1:dim(site_nets)[3], function(i){
	sapply(1:dim(site_nets)[3], function(j){
		compare_networks(site_nets[,,i], site_nets[,,j])
	}, simplify='array')
},simplify='array')

sapply(1:dim(site_nets)[3], function(i){
	compare_networks(full_net, site_nets[,,i])
})




