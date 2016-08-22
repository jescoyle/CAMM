#' Calculate species abundances
#'
#' Calculate species abundances from a community matrix
#'
#' Calculates the abundance of associations or host or symbiont
#' species at each site. The function only counts species
#' present in associations established at microsites and does
#' not count species that are only present in a species' pools.
#'
#' @param comm (required) site x microsite matrix with integers representing  
#' 	the associations present in each microsite
#' @param topo_names (required) host x symbiont association matrix with 
#' 	integers identifying each association (as returned by \code{name_topo})
#' @param type (required) character string indicating whether the function
#' 	should do calculations on host species (\code{'a'}), symbiont
#' 	species (\code{'b'}), or associations (\code{'species'})
#' @param pool site x species matrix for the species pool present at each
#' site (as returned by \code{\link{calc_pool}}). Must be supplied if 
#' calculating for host or symbiont species.
#' @return site x species matrix giving the abundace of each species at each
#' site
#'
#' @export

calc_sa = function(comm, topo_names, type='species', pool=NULL){
	if(type=='species'){
		counts = apply(comm, 1, function(x) table(factor(x, 1:max(topo_names))))
	}
	
	if(type=='a'){
		counts = apply(comm, 1, function(x){
			table(factor(sapply(x, function(xi)	ifelse(xi==0, 0, which(topo_names==xi, arr.ind=T)[1])), 1:nrow(topo_names)))
		})
	}

	if(type=='b'){
		counts = apply(comm, 1, function(x){
			table(factor(sapply(x, function(xi) ifelse(xi==0, 0, which(topo_names==xi, arr.ind=T)[2])), 1:ncol(topo_names)))
		})
	}

	# Do not count species in associations if the species is not present in the species pool
	# This is relevant under facultative mutualism when hosts can establish without symbionts
	if(!is.null(pool)) counts = counts*t(pool)

	t(counts)
}
