#' Compare association networks
#'
#' Compares the assocations shared between two association matrices
#' 
#' The function calculates the number of associations shared between two 
#' association matrices and the number of associations unique to each.
#'
#' @param netA (required) binary species association matrix
#' @param netB (required) binary species association matrix
#' @return a named vector with the number of shared associations 
#' 	(\code{'both'}) and the number of associations unique to each 
#' 	network (\code{'A'} and \code{'B'})
#' 
#' @export

compare_networks = function(netA, netB){

	# Calculate difference between networks
	net_diff = netA - netB

	# Calculate number of shared and unshared associations
	BnotA = sum(net_diff == -1)
	AnotB = sum(net_diff == 1)
	AandB = sum(netA) - AnotB
	
	# Return values
	c(both = AandB, A=AnotB, B=BnotA)
}
