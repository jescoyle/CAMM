#' Calculate abundance
#'
#' Calculates the number of occupied microsites at each site or across the 
#' entire metacommunity
#'
#' @param comm (required) site x microsite matrix with integers representing  
#' 	the associations present in each microsite
#' @return vector of number of occupied microsites at each site or the total
#' 	number of microsites occupied in the entire metacommunity 

#' @describeIn calc_occ Calculate abundance at each site
#' @export
calc_occ = function(comm){
	apply(comm, 1, function(x) sum(x>0))
}

#' @describeIn calc_occ Calculate total abundance for the metacommunity
#' @export
calc_tot_occ = function(comm){
	sum(comm > 0)
}
