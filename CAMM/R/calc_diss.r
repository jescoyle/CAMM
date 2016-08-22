#' Calculate community dissimilarity
#' 
#' Calculate dissimilarity in species composition among sites
#'
#' This function uses \code{\link[vegan]{vegdist}} in the vegan package to
#' calculate community dissimilarity among sites. Dissimilarity can be
#' calculated for associations or for host or symbiont species present
#' in established associations. Any metric implemented by 
#' \code{\link[vegan]{vegdist}} can be used and dissimilarities can be 
#' calculated based on abundances or based on presence/absence 
#' (\code{binary=TRUE}).
#'
#' @inheritParams calc_sa
#' @param metric (required) dissimilarity metric used for calculation
#' 	(see \code{\link{vegdist}}
#' @param binary logical indicating whether dissimilarties should be calculated
#' 	from species abundances (\code{FALSE}, default) or presence/absence
#' 	(\code{TRUE}) 
#' @return matrix of site dissimilarities, as returned by 
#' 	\code{\link[vegan]{vegdist}}
#'
#' @importFrom vegan vegdist
#' @export

calc_diss = function(comm, topo_names, type, metric, binary=F, pool=NULL){
	# Calculate species abundances
	comm_mat = calc_sa(comm, topo_names, type, pool)
	
	# Convert to presence/absence
	if(binary) comm_mat = comm_mat > 0
	
	# Calculate dissimilarities
	vegan::vegdist(comm_mat, method=metric)
}