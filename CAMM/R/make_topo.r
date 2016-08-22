#' Creation of an association network
#'
#' Creates an association network between a given number of host and symbiont
#' species with a defined network topology. 
#' 
#' The function generates a random bipartite network between host and symbiont
#' species which is returned as a matrix of 0s and 1s (indicating whether two 
#' species can associate with one another). Users can specify the network topology
#' as:
#' \describe{
#'	\item{\code{'one2one'}}{each host associates with one symbiont and each symbiont
#'		associates with one host}
#'	\item{\code{'one2many'}}{each host associates with one symbiont and each symbiont
#'		associates with one or more hosts}
#'	\item{\code{'many2many'}}{each host associates with one or more symbionts and
#'		each symbiont associates with one or more hosts}
#' }
#'
#' The specified topology must be compatible with the number of species and 
#' number of links (\code{N_L}) specified. E.g. for \code{'one2one'}, the 
#' number of of host and symbiont species and links must be equal, whereas, for
#' \code{'one2many'} the number of links must be greater than the number of symbionts.
#'
#' @param S_a (required) number of host species
#' @param S_b (required) number of symbiont species
#' @param N_L (required) number of edges (links) in the network
#' @param topology (required) form of network topology. Can be 
#' 	'one2one', 'one2many', 'many2many' (see details).
#' @return a binary matrix with dimension \code{S_a x S_b} indicating links
#' 	between host and symbiont species
#'
#' @export

make_topo = function(S_a, S_b, N_L, topology){
	# Catch errors
	if(N_L < S_a | N_L < S_b) stop('Cannot find partners for all mutualists when N_L < S_a or S_b')
	if(topology=='one2one' & S_a!=S_b) stop('Cannot form one-to-one topology when S_a != S_b')
	if(topology=='one2one' & S_a!=N_L) stop('Cannot form one-to-one topology when S_a != N_L')
	if(topology=='one2many' & N_L <= S_b) stop('Cannot form one-to-many topology when N_L <= S_b')
	if(topology!='many2many' & N_L > S_a) stop('Topology must be many-to-many if N_L > S_a')
	if(N_L > S_a*S_b) stop(paste('Cannot form', N_L, 'associations with', S_a,'and', S_b, 'partners'))

	# Make empty association matrix
	topo = matrix(0, nrow=S_a, ncol=S_b)	

	# Choose one partner for each mutualist 
	partners_a = sample(S_a, S_a, replace=F)
	partners_b = sample(S_b, S_b, replace=F)

	# If topology is many2many add extra links for mutualist a
	if(topology=='many2many') partners_a = c(partners_a, sample(S_a, N_L-S_a, replace=T))

	# If topology is one2many or many2many add extra linkes for mutualist b
	if(topology %in% c('one2many','many2many')) partners_b = c(partners_b, sample(S_b, N_L-S_b, replace=T))

	# Shuffle order
	partners_a = sample(partners_a, replace=F)
	partners_b = sample(partners_b, replace=F)
	

	# Add links to network
	for(i in 1:N_L){
		topo[partners_a[i],partners_b[i]] = 1
	}

	# If there are duplicated links add more random links until at N_L.
	add_links = N_L - sum(topo)
	while(add_links > 0){
		add_a = sample(S_a, add_links, replace=T)
		add_b = sample(S_b, add_links, replace=T)
		for(i in 1:add_links){
			topo[add_a[i], add_b[i]] = 1
		}
		add_links = N_L - sum(topo)
	}

	topo
}