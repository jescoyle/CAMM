#' Calculate species richness
#' 
#' Calculates species richness at each site or across the entire metacommunity
#'
#' Calculates the number of host (\code{type = 'a'}) or symbiont 
#' (\code{type = 'b'}) species present in established associations at each site
#' or the number of associations (\code{type = 'species'}) at each site. 
#'
#' @inheritParams calc_sa
#' @return vector of species richness at each site or the total number of
#' 	species in the entire metacommunity

#' @describeIn calc_rich Calculate site richness
#' @export
calc_rich = function(comm, topo_names, type, pool=NULL){
	sa = calc_sa(comm, topo_names, type, pool)
	apply(sa, 1, function(x) sum(x>0))
}

#' @describeIn calc_rich Calculate total richness
#' @export
calc_tot_rich = function(comm, topo_names, type, pool=NULL){
	sa = calc_sa(comm, topo_names, type, pool)
	sum(colSums(sa) > 0)
}