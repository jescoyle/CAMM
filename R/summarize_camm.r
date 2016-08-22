#' Summarize simulation results
#'
#' Summarizes multiple independent simulation runs
#'
#' This function summarizes a given community descriptor across multiple
#' independent simulations run on the same set of parameters. It extracts
#' the requested descriptor for each simulation run in \code{results}
#' and then calculates the mean, variance, median and 2.5th and 97.5th 
#' quantiles across simulation runs. 
#'
#' The following community descriptors can be requested 
#'(using the paramter \code{what}):
#' \describe{
#'	\item{'S'}{mean site species richness, total species richness, and 
#'		beta diversity}
#'	\item{'N'}{mean site abundance and total abundance}
#'	\item{'Cor_ab'}{correlation between host and symbiont richness}
#'	\item{'cor'}{community ~ environment correlations}
#' }
#' The user must also indicate which community should be summarized
#' useing \code{type}: established associations (\code{'species'}), 
#' or host (\code{'a'}) or symbiont (\code{'b'}) species involved in
#' established associations.
#' 
#' @param results (required) list of simulations results, as returned by
#' 	\code{link{run_sim_N}}
#' @param what (required) character string indicating which community 
#' 	descriptor to summarize (see details)
#' @param type character string indicating which community to summarize
#' 	(see details: required if \code{what} is \code{'S'}, \code{'N'}, or 
#' 	\code{'cor'})
#' @return an array
#'
#' @export

summarize_camm = function(results, what, type=NA){
	# Determine whether results have multiple time points
	tp = length(dim(results[[1]][[1]]))==3

	if(tp){
		# Richness in each community
		if(what=='S'){
			get_vars = paste(c('S','Stot','Beta'),type,sep='_')
			richness = sapply(results, function(x) x[[1]][,get_vars,], simplify='array')
			return_stats = apply(richness, 1:3, function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response','time')
		}
		
		# Abundance
		if(what=='N'){
			abun = sapply(results, function(x) x[[1]][,c('N', 'Ntot'),], simplify='array')
			return_stats = apply(abun, 1:3, function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response','time')
		}
	
		# Correlation between host and symbiont richness
		if(what=='Cor_ab'){
			abun = sapply(results, function(x) x[[1]][,c('Cor_ab'),], simplify='array')
			return_stats = apply(abun, 1:2, function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','time')
		}
		
		# Community ~ environment correlations
		if(what=='cor'){
			corr_mat = sapply(results, function(x) x[[2]][,type,,,,], simplify='array')
			maxD = length(dim(corr_mat))
			return_stats = apply(corr_mat, 1:(maxD-1), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','env','measure','time')
		}
	} else {
		# Richness in each community
		if(what=='S'){
			get_vars = paste(c('S','Stot','Beta'),type,sep='_')
			richness = sapply(results, function(x) x[[1]][,get_vars], simplify='array')
			return_stats = apply(richness, c(1,2), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response')
		}

		# Abundance
		if(what=='N'){
			abun = sapply(results, function(x) x[[1]][,c('N', 'Ntot')], simplify='array')
			return_stats = apply(abun, c(1,2), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','response')
		}

		# Correlation between host and symbiont richness
		if(what=='Cor_ab'){
			corab = sapply(results, function(x) x[[1]][,c('Cor_ab')], simplify='array')
			return_stats = apply(corab, 1, function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary')
		}

		# Community ~ environment correlations
		if(what=='cor'){
			corr_mat = sapply(results, function(x) x[[2]][,type,,,], simplify='array')
			maxD = length(dim(corr_mat))
			return_stats = apply(corr_mat, 1:(maxD-1), function(x) c(mean=mean(x, na.rm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975), na.rm=T)))
			names(dimnames(return_stats)) = c('stat','summary','env','measure')
		}
	}

	return_stats
}
