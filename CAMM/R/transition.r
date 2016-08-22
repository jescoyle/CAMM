#' Transition microsite
#'
#' Stochastically transitions a microsite in a community based on transition
#' probabilities.
#'
#' @param probs (required) named vector of probabilities that a microsite will
#' 	take on a given value in the next timestep. Names are the state values.
#' @return an integer representing the new state

transition = function(probs){
	as.numeric(sample(names(probs), 1, prob=probs))
}