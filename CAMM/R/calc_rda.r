#' Community-environment RDA
#'
#' Calculate community ~ environment RDA
#' 
#' The function calculates the correspondence between community composition
#' and a given environmental data matrix (\code{env}) for established
#' associations or host or symbiont species in established associations
#' by performing a constrained ordination (\code{\link[vegan]{rda}}) on 
#' Hellinger-transformed species abundances (as implemented 
#' by \code{\link[vegan]{decostand}}).
#' The function returns the adjusted eqn{R^2} as a measured of correspondence
#' between community structure and environment.
#' 
#' @inheritParams calc_sa
#' @param env (required) matrix of environmental variables across sites
#' @param binary logical indicating whether rda should be performed on
#' 	species abundances (\code{FALSE}, default) or presence/absence
#' 	(\code{TRUE})
#'
#' @seealso \code{\link[vegan]{rda}} for constrained ordination
#'
#' @importFrom vegan decostand rda RsquareAdj

calc_rda = function(comm, topo_names, type, env, binary=F, pool=NULL){
	# Calculate community matrix of species abudances
	comm_mat = calc_sa(comm, topo_names, type, pool)

	# Convert to presence/absence if indicated
	if(binary) comm_mat = comm_mat > 0
	
	# Transform abundances
	comm_std = vegan::decostand(comm_mat, 'hellinger')

	# Conduct constrained ordination
	ord = vegan::rda(comm_std, env)

	# Return the R-square
	#RsquareAdj(ord)$adj.r.squared
	vegan::RsquareAdj(ord)$r.squared
}
