#' Calculate association network
#'
#' Calculates associations network from observed established communities
#' 
#' This function calculates the observed species association network from 
#' a metacommunity (\code{comm}). It can return either a single matrix
#' across the entire metacommunity or an array of matrices, one for each 
#' site, if \code{bysite=TRUE}.
#'
#' @param comm (required) site x microsite matrix with integers representing  
#' 	the associations present in each microsite
#' @param topo_names (required) host x symbiont association matrix with 
#' 	integers identifying each association (as returned by \code{name_topo})
#' @param bysite logical indicating whether networks should be calculated for
#' 	each site in \code{comm}
#' @return a binary matrix or array representing the association network(s)
#'
#' @export



calc_network = function(comm, topo_names, bysite=F){

	# Calculate presence of associations as each site
	abuns = calc_sa(comm, topo_names, type='species')>0
	
	# Calculate network for each site
	topos = sapply(1:nrow(abuns), function(i){
		x = abuns[i,]
		this_topo = 1*(topo_names>0)
		this_topo[which(!x)] = 0
		this_topo
	}, simplify='array')

	# Sum across sites, if requested
	if(!bysite) topos = apply(topos, 1:2, function(x) 1*(sum(x)>0))

	topos
}