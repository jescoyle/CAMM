# CAMM
The Community Assembly of Mutualists Model (CAMM) is a general simulation model of community assembly for multiple host and symbiont species. The purpose of building this model was to explore the consequences of obligate mutualism on diversity patterns along environmental gradients. This model extends previous niche-based community assembly models to two communities (hosts and symbionts) linked by an association network that defines which species pairs facilitate each others’ establishment. CAMM is a spatially implicit patch-occupancy model that represents an open meta-community with patches jointly distributed along two environmental gradients. Patches are comprised of a fixed number of microsites, each of which may be occupied by a host and its symbiont. Each patch also has a pool of unassociated host and symbiont species. The state of each microsite follows a discrete-time Markov process with transition probabilities determined by the presence of species in the pool of established an unassociated mutualists. The probability of species occuring in a patch's pool is determined by the match between species' niches and the environmental conditions in the patch as well as species' global abundances. Thus, the primary process built into the model is environmental filtering. Competition only occurs as microsite pre-emption and there are no differences in mortality among species or patches.

The simulation currently operates as follows:

1. Randomly initialize pools of host and symbiont species with 2-dimensional environmental niches, a landscape of patches with envionmental conditions, and an association network between hosts and symbiont. The simulation starts from an empty metacommunity with N[C] patches each with space for up to N individuals.
2. During each time step the following operations occur in each patch in the following order:
  1. Unassociated hosts and symbionts die with fixed probability.
  2. Host and symbiont propagule disperse (independently) into patches probabilistically according to fixed global species distributions.
  3. A transistion probability matrix is calculated that gives the probability that an empty microsite or established mutualism will stay the same or change status in this timestep. In the current simulation, the probability that an empty space is colonized by a host and symbiont is a function of the presence of hosts and symbionts in the patch's pool of species, the degree to which hosts depend on symbionts for establishment and the probability that a host will associate with a given symbiont (currently the same for all associations in all environments). The probability that colonized microsite will transition to an empty microsite is determined by a fixed death rate and the probability that colonized microsites will transition directly to a different set of mutualists is set to 0 (no competition).
  4. Each microsite stochastically changes to a new state according to the transition matrix.
3. The simulation is run for a fixed number of timesteps to reach an equllibrium.
4. Descriptive statistics are calculated to summarize patterns of commmunity variation along the two environmental gradients.
5. Steps 2-4 is repeated several times on the same initial metacommunity ("multiple chains") and descriptive statistics are averaged.
6. Step 5 is repeated multiple times to generate a distribution of community descriptive statistics for a given set of simulation parameters.


## Files
### ./R : Contains scripts for running the simulation in R
+ parameter_file.R : A file containing parameter values for one simulation run.
+ simulation_functions.R : Contains functions for initializing and running the simulation.
+ run_camm.R : The script used to run one implementation of the simulation.
+ analysis_function.R : Contains functions for analyzing simulation results.
