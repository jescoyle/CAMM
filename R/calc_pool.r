#' Calculate species pool
#'
#' Determines which species are present in a community.
#'
#' @param comm (required) metacommunity matrix indicating which 
#' 	association is present in microsite of each site. Rows are sites
#' 	and columns are microsites/individuals.
#' @param topo_names (required) matrix representing association network 
#' 	where integers label each association, as produced by 
#' 	\code{\link{name_topo}}. Labels should match those used in \code{comm}.
#' @param partner (required) integer indicating whether the function should
#' 	calculate the species pool for hosts (\code{1}) or symbionts (\code{2}).
#' @return site x species matrix indicating species presence in each site

calc_pool = function(comm, topo_names, partner){
	# Find number of communities and species
	N_C = nrow(comm)
	N_S = dim(topo_names)[partner]

	# For each row of comm (representing a community) find the associations present and find the partners present
	pool = matrix(0, nrow=N_C, ncol=N_S)
	for(i in 1:N_C){

		# Find the associations present
		assocs = unique(comm[i,])
		assocs = assocs[assocs!=0]

		# Find the partners present
		for(x in assocs){
			j = which(topo_names==x, arr.ind=T)[,partner]

			# Assign 1 for presence in community i
			pool[i,j] = 1
		}
	}

	# Return partner pool
	pool
}