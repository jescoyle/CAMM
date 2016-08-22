#' Calculate transition probability matrices
#'
#' Calculates the probabilities that microsites will change states
#' for each site in a metacommunity.
#' 
#' For each site in a metacommunity, the function calculates a transition
#' matrix containing the probability that a microsite will changes states
#' during a single timestep. Microsites can be empty (0) or contain an
#' established host and symbiont (an 'association', given as an integer, 
#' \code{i}). For each potential association (given in \code{topo_names}), 
#' the function determines whether both host and symbiont are present in 
#' the species pools (\code{poolA} and \code{poolB}). If the host is not 
#' present, the association cannot establish (\code{P(0 -> i) = 0}). If both 
#' mutualists are present then the probability of an association forming is
#' given by \code{assoc_probs}, whereas if the symbiont is not present then
#' this probability is multiplied by \code{1 - omega}, which allows for 
#' facultative mutualism. Then, the probability that no association establishes
#' (\code{P(0 -> 0)}) is calculated as the product of the probabilities of  
#' each potential association failing to form. The remaining probability
#' is divided proportionally among each of the previously calculated
#' association probabilities (since the rows of a transition matrix must sum
#' to \code{1}). The probability that an association dies (\code{P(i -> 0)})
#' is a constant that is the same for all associations (\code{mort_rate}) 
#' and there is no displacement of one association by another
#' (\code{P(i -> j) = 0}). However, the function could be modified to allow for
#' competitive displacement.
#'
#' @param sites (required) a 2 column matrix giving environmental conditions
#' 	at each site
#' @param niches_a (required) array of niche optima and breadths for host 
#' 	species on each environmental axis (\code{[species, optima/breadth,
#' 	environment]}) as returned by \code{\link{make_niches_gsad}}.
#' @param niches_b (required) array of niche optima and breadths for symbiont 
#' 	species on each environmental axis (\code{[species, optima/breadth,
#' 	environment]}) as returned by \code{\link{make_niches_gsad}}.
#' @param topo_names (required) matrix representing association network 
#' 	where integers label each association, as produced by 
#' 	\code{\link{name_topo}}.
#' @param poolA (required) binary site x species matrix indicating which host 
#' 	species are present in each site's species pool, but which are not 
#' 	necessarily established.
#' @param poolB (required) binary site x species matrix indicating which symbiont 
#' 	species are present in each site's species pool, but which are not 
#' 	necessarily established.
#' @param mort_rate (required) probability that an association dies
#' @param assoc_probs (required) an site x association matrix giving the 
#' 	probability that an association forms at each site if both partners 
#' 	are present. Has dimensions \code{[sites, associations]}.
#' @param omega (required) number controling the strength of mutualism.
#' 	\code{1-omega} is the probability that a host will establish without 
#' 	its symbiont.
#' @return array of transition matrices, one for each site. 
#' 	Has dimension \code{[association, association, site]}.

calc_probs = function(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs, omega){
	# Calculate number of associations and communities
	S_a = nrow(topo_names)
	S_b = ncol(topo_names)
	N_L = max(topo_names)
	N_C = nrow(sites)

	# Create empty array that holds transition matrices for each site
	Tarr = array(NA, dim=c(N_L+1, N_L+1, N_C), dimnames=list(0:N_L, 0:N_L, 1:N_C))

	# A 2 x N_C x N_L matrix indicating whether partner A or B is present in each community
	AB_present = sapply(1:N_L, function(x){
		partners = which(topo_names==x, arr.ind=T)
		rbind(poolA[,partners[1]],poolB[,partners[2]])
	}, simplify='array')	


	# Calculate probabilities of establishment 0 -> i
	establish = apply(AB_present, c(2,3), function(x){
		# If partner A is absent P = 0
		# Otherwise, if partner B is absent, P = 1-omega
		# Or, if partner B is also present, P = 1
		ifelse(x[1]==0, 0, ifelse(x[2]==0, 1-omega, 1))
	})*assoc_probs 
	
	# Calculate probability that no species establishes
	no_establish = apply(establish, 1, function(p) prod(1-p))

	# Scale each P(0 -> i) so that rows sum to 1. A species has P=establish only if it is the only species that can establish. 
	establish = establish/ifelse(rowSums(establish)>0, rowSums(establish),1)*(1-no_establish)# scale establishment probabilities so that 0 
	Tarr[1,,] = t(cbind(no_establish, establish))
	
	# Add probability of death
	# For now, constant for any established species
	Tarr[2:(N_L+1),1,] = mort_rate	
	
	# For each species add transition probabilities
	for(i in 1:N_L){
		# Probability of competitive displacement
		# For now, no competitive displacement
		# With competition, entries should be the product of:
			# probability that the competitor can establish 
			# a (environmentally dependent) competitive  ability coefficient
			# scaled to sum=1 if sum > 1
		compete = matrix(0, N_C, N_L-1)
		colnames(compete) = (1:N_L)[1:N_L != i]
		compete = (1-Tarr[i+1,1,])*compete # Multiply competitive probabilities by P(i survives)

		# Probability of remaining
		remain = 1 - Tarr[i+1,1,]- rowSums(compete)

		Tarr[i+1,i+1,] = remain
		Tarr[i+1, colnames(compete),] = t(compete)		

	}
	# Check that rows sum to 1
	if(sum(round(apply(Tarr, c(1,3), sum),10)!=1)>0) stop('Transition matrix rows do not sum to 1')

	# Return transition array
	Tarr
}
