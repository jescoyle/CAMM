#' Death
#'
#' Random mortality of species in communities that are not
#' established in a microsite as a mutualistic association.
#'
#' @param comm (required) metacommunity matrix indicating which 
#' 	association is present in microsite of each site. Rows are sites
#' 	and columns are microsites/individuals.
#' @param topo_names (required) matrix representing association network 
#' 	where integers label each association, as produced by 
#' 	\code{\link{name_topo}}. Labels should match those used in \code{comm}.
#' @param pool (required) binary site x species matrix indicating which species
#' 	 are present in each site's species pool, but which are not necessarily
#' 	established.
#' @param mortality (required) probability that an unassociated mutalist dies
#' @param partner (required) integer indicating whether the function
#' 	applies to hosts (\code{1}) or symbionts (\code{2}).
#' @return binary site x species matrix indicating which species remain
#' 	in the species' pools at each site


# A function that causes random mortality of mutualists in a community. Mutualists involved in an association cannot die.
#	comm = matrix indicating which association is present in each position of a community
#	topo_names = bipartite mutualist network with integers labeling associations
#	pool = site X mutualist binary matrix indicating which mutualist species are present
# 	mortality = probability that an unassociated mutualist dies
#	partner = integer indicated which partner this applies to
die = function(comm, topo_names, pool, mortality, partner){
	
	# Find pool of mutualists that are not in an association and can die
	comm_pool = calc_pool(comm, topo_names, partner)
	not_associated = which(pool!=comm_pool) # also finds species in associations that are not in the pool (when omega < 1)
	
	# For each unassociated partner in each community do random mortality
	survive = runif(length(not_associated)) > mortality
	pool[not_associated] = pool[not_associated]*survive # species in associations that are not in the pool remain not in pool while species in pool but not associations may die

	# Return new pool
	pool
}