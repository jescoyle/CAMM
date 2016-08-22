#' Community-environment correlations
#'
#' Calculates multiple community-environment correlations
#'
#' This summary function calculates the correlation between each
#' environmental axis (provided in \code{env}) and species richness, abundance
#' and community structure for established associations and host and symbiont
#' species involved in established associations. For each environmental axis and 
#' dissimilarity metric given in \code{metric} it calculates the correlation 
#' between a community dissimilarity matrix (see \code{link{calc_diss}}) and
#' the environmental distance matrix. Finally, the function conducts a 
#' constrained ordination along each environmental axis (see  
#' \code{link{calc_rda}}) and returns the adjusted \eqn{R^2}. These statistics
#' are returned as a 4-dimensional array where the first dimension indicates 
#' the community examined (associations: \code{'species'}, hosts: \code{'a'}, 
#' symbionts: \code{'b'}), the second dimension indicates which environmental
#' axis was used, the third dimension indicates the statistic (richness:
#' \code{'S'}, abundance: \code{'N'}, RDA: \code{'rda'}, or correlation with
#' community dissimilarity metric), and the fourth dimension indicates whether
#' the calculations are based on abundances (\code{FALSE} or presence/absence
#' (\code{TRUE}). 
#'
#' @param comm (required) site x microsite matrix with integers representing the 
#' 	associations present in each microsite
#' @param topo_names (required) host x symbiont association matrix with integers
#' 	identifying each association (as returned by \code{name_topo})
#' @param env (required) site x env matrix with environmental conditions at 
#' 	each site
#' @param metric vector of dissimilarity metric (see 
#' 	\code{\link[vegan]{vegdist}} for options). Defaults to \code{NULL} for none.
#' @param binary lvector of logicals indicating whether calculations should be
#' 	performed on species abundances (\code{FALSE}, default) or presence/absence
#' 	(\code{TRUE}) or both
#' @param pools (required) named list of pools of hosts as symbionts present at
#' 	each site. Pools are generated from the community matrix by 
#' 	\code{calc_pool}. The host pool should be named \code{'a'} and the symbiont 
#' 	pool named \code{'b'}. 
#' @return array with dimensions \code{[community type, environmental axis, 
#' 	measurement, abundances vs presence/absence]} (see details)
#'
#' @export

calc_envcorr = function(comm, topo_names, env, metric=NULL, binary, pools=NULL){
	# Count the number of environmental axes
	Nenv = ifelse(is.null(dim(env)), 1, ncol(env))

	# Define array to hold statistics
	corr_mat = array(NA, dim=c(3, Nenv, length(metric)+3, length(binary)), 
		dimnames=list(community=c('species','a','b'), env=1:Nenv, measure=c('S','N','rda', metric), binary=binary))

	for(k in binary){
	for(type in c('species','a','b')){
	for(i in 1:Nenv){
		# Get environmental axis
		if(Nenv==1){ X = env } else { X = env[,i]}
		
		# Calculate correlation with species richness
		# NOTE: DOES NOT IMPLEMENT ABUNDANCE-WEIGHTED SPECIES DIVERSITY, SAME CALCULATION FOR BINARY=T/F
		if(type=='species'){
			corr_mat[type, i, 'S', as.character(k)] = cor(X, calc_rich(comm, topo_names, type))
		} else {
			corr_mat[type, i, 'S', as.character(k)] = cor(X, calc_rich(comm, topo_names, type, pools[[type]]))
		}

		# Calculate correlation with abundance
		# NOTE: THIS WILL BE THE SAME REGARDLESS OF TYPE OR BINARY
		corr_mat[type, i, 'N', as.character(k)] = cor(X, calc_occ(comm))

		# Calculate RDA
		if(type=='species'){
			corr_mat[type, i, 'rda', as.character(k)] = calc_rda(comm, topo_names, type, X, k)
		} else {
			corr_mat[type, i, 'rda', as.character(k)] = calc_rda(comm, topo_names, type, X, k, pools[[type]])
		}
		
		# Calculate correlation for each dissimilarity metric
		if(length(metric)>0){
			for(j in metric){
				if(type=='species'){
					comm_diss = calc_diss(comm, topo_names, type, j, k)
				} else {
					comm_diss = calc_diss(comm, topo_names, type, j, k, pools[[type]])
				}
				env_dist = dist(X)
				corr_mat[type, i, j, as.character(k)] = cor(comm_diss, env_dist, use='pairwise.complete.obs')
			}
		}

	}}}
	
	corr_mat
}
