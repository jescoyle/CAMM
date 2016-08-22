#' Identifies links in an association network
#'
#' Labels each link in an association network with a unique integer. 
#' 
#' @param topo (required) association network supplied as a binary matrix
#' @return matrix with the same structure as \code{topo} in which all 1s
#' have been replaced by integers


name_topo = function(topo){
	topo_labeled = t(topo)
	N_L = sum(topo_labeled==1)
	topo_labeled[topo_labeled==1] = 1:N_L

	t(topo_labeled)
}
