# CAMM
The Community Assembly of Mutualists Model (CAMM) is a general simulation model of community assembly for multiple host and symbiont species. The purpose of building this model was to explore the consequences of obligate mutualism on diversity patterns along environmental gradients. This model extends previous niche-based community assembly models to two communities (hosts and symbionts) linked by an association network that defines which species pairs facilitate each others’ establishment. CAMM is a spatially implicit patch-occupancy model that represents an open meta-community with patches jointly distributed along two environmental gradients. Patches are comprised of a fixed number of microsites, each of which may be occupied by a host and its symbiont. Each patch also has a pool of unassociated host and symbiont species. The primary process built into the model is environmental filtering. Competition only occurs as microsite pre-emption and there are no differences in dispersal among species or patches.

The simulation currently operates as follows:
1. Randomly initialize pools of host and symbiont species with 2-dimensional environmental niches, a landscape of patches with envionmental conditions, and an association network between hosts and symbiont.
2. 



## Files
### ./R : Contains scripts for running the simulation in R
+ parameter_file.R : A file containing parameter values for one simulation run.
+ simulation_functions.R : Contains functions for initializing and running the simulation.
+ run_camm.R : The script used to run one implementation of the simulation.
+ analysis_function.R : Contains functions for analyzing simulation results.
