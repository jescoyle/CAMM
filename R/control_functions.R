## This script holds functions to control the community assembly of mutualists model


# A function that starts multiple runs of a simulation
# parm_file = a file containing parameter values
# nruns = number of simulation runs to start
# nparallel = if run in parallel, number of cores
# sim_parms = named list of parameter controling when simulations should stop
#	components are: sim_mode and reps or nchains (depending on sim_mode)
# simID = string that identifies this set of simulations runs
# save_start = flag indicating whether each initial set of niches and topologies should be saved
run_camm_N = function(parm_file, nruns, nparallel=1, sim_parms, simID, save_start=F){
	
	foreach
		runID = paste0(simID, i)
		metacomm = initialize_camm(parm_file, save_start, runID=runID)
		sim_results = run_camm(metacomm, sim_mode=sim_parms$sim_mode, reps=sim_parms$reps, nchains=sim_parms$nchains)
		


} 



