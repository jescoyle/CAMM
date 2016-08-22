#' Creation of sites 
#'
#' Randomly generates sites with two potentially correlated environmental
#' conditions. Conditions are drawn from a multivariate normal distribution
#' with mean = 0 and variance = 1. 
#'
#' @param N_C (required) number of sites
#' @param rho_z correlation between two environmental conditions
#' @return a matrix of values with dimension \code{N_C x 2}
#'
#' @export
#' @importFrom MASS mvrnorm

make_sites = function(N_C, rho_z=0){
	sd = c(1,1)
	mu = c(0,0)

	# Define covariance matrix
	S = matrix(c(sd[1],rho_z,rho_z,sd[2]), nrow=2, ncol=2)

	# Generate correlated gaussian rv using MASS package
	sites = MASS::mvrnorm(n=N_C, mu=mu, Sigma=S)
	
	# Return sites
	sites
}
