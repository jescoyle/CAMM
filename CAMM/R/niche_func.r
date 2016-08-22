#' Calculate survival probability
#' 
#' Calculates the probability that a propagule will establish
#' under specific environmental conditions
#'
#' The probability that a propagule will establish is the product
#' of the probabilities that the propagule will establish under each
#' environmental condition. (E.g. it assumes independent effects of the
#' environmental axes.) Probabilities are calculated using a Gaussian curve
#' whose mean is the niche optimum and whose standard deviation is the 
#' niche breadth. Curves are not standardized to be pdfs; thus, each curve
#' reaches a maximum of 1 when the environment exactly matches a species 
#' optimum.
#' 
#' @param x (required) vector of environmental conditions
#' @param mu (required) vector of niche optima (matching x)
#' @param sigma (required) vector of niche breadths (matching x)
#' @return the probability of establishment

niche_func = function(x, mu, sigma){
	N = length(x)
	pvals = sapply(1:N, function(i) 1/exp((x[i]-mu[i])^2/(2*sigma[i]^2)))

	prod(pvals)
}
