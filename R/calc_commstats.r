#' Calculate community statistics
#'
#' Calculate diversity and abundance of hosts, symbionts and associations
#'
#' This function calculates abundance and diversity of hosts, symbionts
#' and associations that are established at sites. If a species is present
#' in a site's species pool, but not in an association at a microsite it 
#' will not be counted. Three diversity statistics are calculated for hosts
#' (\code{'a'}), symbionts (\code{'b'}) and associations (\code{'species'}):
#' \itemize{
#'	\item \strong{\code{S}} : species richness (number of species present
#'		at each site), see \code{\link{calc_rich}}
#'	\item \strong{\code{Stot}} : total species richness across all sites,
#'		see \code{\link{calc_tot_rich}}
#'	\item \strong{\code{Beta}} : beta diversity (total richness divided
#'		by mean richness at each site)
#' }
#' The function also returns the number of microsites occupied at each site
#' (\code{N}), the total number of microsites occupied across the entire
#' metacommunity (\code{Ntot}), and the pearson correlation between 
#' species richness of hosts and symbionts across sites (\code{Core_ab}).
#' Community statistics are returned as a dataframe with each row representing
#' a site.
#'
#' @param comm (required) site x microsite matrix with integers representing the 
#' 	associations present in each microsite
#' @param topo_names (required) host x symbiont association matrix with integers
#' 	identifying each association (as returned by \code{name_topo})
#' @param pools (required) named list of pools of hosts as symbionts present at
#' 	each site. Pools are generated from the community matrix by 
#' 	\code{calc_pool}. The host pool should be named \code{'a'} and the symbiont 
#' 	pool named \code{'b'}.
#' @return a dataframe of community statistics (see details) with rows 
#' 	representing sites
#' 
#' @export

calc_commstats = function(comm, topo_names, pools){

	# Calculate richness
	S_species = calc_rich(comm, topo_names, 'species')
	S_a = calc_rich(comm, topo_names, 'a', pools[['a']])
	S_b = calc_rich(comm, topo_names, 'b', pools[['b']])
	
	# Calculate total richness
	Stot_species = calc_tot_rich(comm, topo_names, 'species')
	Stot_a = calc_tot_rich(comm, topo_names, 'a', pools[['a']])
	Stot_b = calc_tot_rich(comm, topo_names, 'b', pools[['b']])

	# Calculate beta diversity
	Beta_species = Stot_species / mean(S_species, na.rm=T)
	Beta_a = Stot_a / mean(S_a, na.rm=T)
	Beta_b = Stot_b / mean(S_b, na.rm=T)

	# Calculate correlation between richness of a and b
	Cor_ab = cor(S_a, S_b, use='complete.obs')

	# Calculate total abundance
	N = calc_occ(comm)
	Ntot = calc_tot_occ(comm)

	# Return dataframe
	# NOTE: tots will be the same across all communities
	data.frame(S_species, Stot_species, Beta_species, S_a, Stot_a, Beta_a, S_b, Stot_b, Beta_b, N, Ntot, Cor_ab)
}