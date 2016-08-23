## This script defines state variables used to run the community assembly of mutualists model (CAMM)

## Name simulation run
runID = 'base'

## Define mutualist names
a_name = 'myco' # The primary partner or host
b_name = 'photo' # The secondary partner or symbiont

## Define number of taxa for mutualists a and b
S_a = 10
S_b = 3

## Define number and size of communities
N_C = 6 # number of sampled communities
N = 10 # number of individuals in a single community

## Parameters controling environmental variation
rho_z = 0 # correlation of two environmental variables across sites

## Parameters controling niches
# Niche optima - drawn from a uniform distribution
mu_a1 = 1 # maximum niche optimum of mutualist a for environmental variable 1
mu_a2 = 1 # maximum niche optimum of mutualist a for environmental variable 2
rho_a = 0 # correlation of niche optima for 2 environmental variables of mutualist a
mu_b1 = 1 # maximum niche optimum of mutualist b for environmental variable 1
mu_b2 = 1 # maximum niche optimum of mutualist b for environmental variable 2
rho_b = 0 # correlation of niche optima for 2 environmental variables of mutualist b

# Niche breadth - drawn from a gamma distribution, alphas may need to be the same and >1 for generating niches 
sigma_a1 = 10 # mean niche breadth of mutualist a for environmental variable 1
sigma_a2 = .5 # mean niche breadth of mutualist a for environmental variable 2
alpha_a1 = 10 # shape parameter of gamma distribution of niche breadths for mutualist a for environmental variable 1
alpha_a2 = 10 # shape parameter of gamma distribution of niche breadths for mutualist a for environmental variable 2
r_a = 0 # covariance of niche breadths of 2 environmental variables for mutualist a: negative -> tradeoffs
sigma_b1 = .5 # mean niche breadth of mutualist b for environmental variable 1
sigma_b2 = 10 # mean niche breadth of mutualist b for environmental variable 2
alpha_b1 = 10 # shape parameter of gamma distribution of niche breadths for mutualist b for environmental variable 1
alpha_b2 = 10 # shape parameter of gamma distribution of niche breadths for mutualist b for environmental variable 2
r_b = 0 # covariance of niche breadths of 2 environmental variables for mutualist b: negative -> tradeoffs

## Parameters controling global abundance distribution (differential dispersal)
# Shape of the SAD

gsad_dist_a = list(type='same', maxN = 2, P_maxN = 0.999, corr=matrix(0, nrow=2, ncol=2))
gsad_dist_b = list(type='same', maxN = 2, P_maxN = 0.999, corr=matrix(0, nrow=2, ncol=2))

# Lists of parameters for each mutualist
nicheparms_a = list(mu = c(mu_a1, mu_a2), rho = rho_a, sigma = c(sigma_a1, sigma_a2), alpha = c(alpha_a1, alpha_a2), r=r_a)
nicheparms_b = list(mu = c(mu_b1, mu_b2), rho = rho_b, sigma = c(sigma_b1, sigma_b2), alpha = c(alpha_b1, alpha_b2), r=r_b)

## Parameters controling mutualistic network

# Topology
N_L = S_a # number of links: must be >= S_a, if > then must be topology = many2many'
topology = 'one2many' # links from mutualist a to mutualist b: 'one2one','one2many','many2many'

# Strength of mutualism
# Range from 0 (no mutualism) to 1 (obligate mutualism)
# Currently affects the probability of establishment when partner A is present and B is not
# P(establish| A & !B) = (1 - omega)*P(association)
omega = 0 

# Environmental dependence of association
# 	a function that calculates the probability of association given environmental conditions at a site and optimal niches of mutualists
# 	each argument is a vector of two numbers, one for each environmental axis
assoc_prob_func = function(env, niche_a, niche_b){
	0.9	
}

# Benefits of association
mort_rate = 0.05 # mortality rate of associations: COULD CHANGE TO VARY AMONG ASSOCIATIONS
mort_rate_a = 10 # mortality rate of unassociated mutualist a relative to association
mort_rate_b = 10 # mortality rate of unassociated mutualist b relative to association













