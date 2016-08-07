#' Creation of metacommunity
#'
#' Create a set of \code{N_C} communities, each with \code{N} microsites
#' or individuals. 
#' 
#' A metacommunity is represented by a \code{N_C x N} matrix where rows give
#' the host-symbiont pair present at each site. Integers in the matrix
#' correspond to host-symbiont pairs named in the association network topology
#' matrix. Currently implements two methods for starting a metacommunity which
#' is specified by \code{fill}:
#' \decribe{
#'	\item{\code{'empty'}}{initial metacommunity has no species present}
#'	\item{\code{'random'}}{initial metacommunity community populated by randomly 
#'		sampling species pairs until it is full}
#'
#' @param N_C (required) number of communities (sites)
#' @param N (required) number of microsites (or individuals) in each community
#' @param N_S number of species pairs
#' @param fill character string indicating how communities should be filled.
#' 	Defaults to \cpde{'empty'}. See details.

make_comm = function(N_C, N, N_S=0, fill='empty'){
	# If starting with empty communities
	if(fill=='empty') comm = matrix(0, nrow=N_C, ncol=N)
	
	# If filling communities with random species pairs
	if(fill=='random') comm = matrix(sample(0:N_S, N_C*N, replace=T), N_C, N)	
	
	comm
}
