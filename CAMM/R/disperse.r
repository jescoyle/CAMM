#' Dispersal
#' 
#' Stochastic dispersal of species into communities based on the 
#' global species abundance distribution
#'
#' The global species abundance distribution \code{gsad} gives
#' the number of propagules of each species arriving at each site
#' in one time step. The probability that a species establishes 
#' is the probability that at least on propagule establishes, given
#' a per-propagule establishment probability calculated by 
#' \code{\link{niche_func}} based on the match between species niches
#' and environmental conditions in the site. The function assumes that 
#' hosts and symbionts disperse independently (no co-dispersal), since it
#' can only be applied to one or the other and not both concurrently.
#'
#' @param sites (required) a 2 column matrix giving environmental conditions
#' 	at each site
#' @param niches (required) array of niche optima and breadths for each species
#' 	and each environmental axis (\code{[species, optima/breadth, environment]})
#' 	as returned by \code{\link{make_niches_gsad}}.
#' @param pool (required) binary site x species matrix indicating which species
#' 	 are present in each site's species pool, but which are not necessarily
#' 	established.
#' @param gsad vector of integers giving the number of propagules for each species
#' 	arriving at a site in a timestep
#' @return binary site x species matrix indicating which species are now
#' 	in the species' pools at each site

disperse = function(sites, niches, pool, gsad){
	
	# Find number of communities and species
	N_C = nrow(pool)
	N_S = ncol(pool)

	# For each species and each site calculate the probability of establishment based on joint gaussian niches on each environmental variable
	probs = matrix(0,N_C, N_S)	

	for(i in 1:N_C){
	for(j in 1:N_S){
		p = niche_func(sites[i,], niches[j,'mu',], niches[j,'sigma',])
		probs[i,j] = 1 - dbinom(0, gsad[j], p) # Probability of at least one propagule establishing with propagule pressure equal to global abundance
	}}

	# Generate random dispersal
	rands = runif(N_C*N_S) 
	emmigrants = rands <= probs
	
	# Add emmigrants to existing pool
	pool = pool + emmigrants > 0
	pool = matrix(as.numeric(pool), N_C, N_S)

	# Return new pool of partners
	pool
}

