#' Creation of species niches and global abundances
#'
#' Concurrently generates random 2-D gaussian niches and a species abundance
#' distribution that is correlated with niche parameters.
#' 
#' This function randomly generates niche optima and breadths for \code{N_S} species
#' along two envronmental axes as well as a global species abundance
#' distribution. Niche optima and breadths may be correlated with global 
#' abundance along one or both environmental axes. In additon, niche optima and 
#' breadths may also be correlated across axes. These correlations are
#' controlled by parameters in the named lists \code{nicheparmes} and 
#' \code{distribution} and are accomplished through Gaussian copulas. This method
#' samples five variables (niche optima on axes 1 and 2, niche breadths on axes
#' 1 and 2, and species global abundance) from a multivariate normal 
#' distribution with covariance matrix defined by \code{r} (correlation between
#' niche optima across environmental axes), \code{rho} (correlation between 
#' niche breadths across environmental axes), and \code{corr} (correlation 
#' between global abundance and each niche parameter, given as a 2x2 matrix 
#' where columns refer to optima and breadths and rows refer to environmental 
#' axes). Then, niche optima are transformed to a uniform distribution on the  
#' interval \code{[-mu, mu]} for each environmental axis. Niche breadths
#' are transformed to a gamma distribution with mean \code{sigma} and shape 
#' parameter \code{alpha} for each environmental axis. Species abundances are 
#' transformed to the distribution given in \code{type}. The following forms
#' of global species abundance distributions are supported and are described 
#' along with thie requried parameters.
#' \describe{
#'	\item{\code{'same'}}{All species have the same abundance.}
#'	\item{\code{'uniform'}}{Uniform distribution on \code{[1, maxN]}.
#'		Can also specify endpoint with \code{p1}.}
#'	\item{\code{'power'}}{Power-law distribution with power \code{p1}.
#'		Can alternatively specify \code{maxN} and \code{P_maxN} to find  
#'		\code{p1} corresponding to a distribution with a probability of  
#'		\code{P_maxN} of exceeding \code{maxN}.}
#'	\item{\code{'logseries'}}{Logseries distribution where \code{p1} is
#'		Fisher's alpha and \code{p2} is N. Can alternatively specify 
#'		\code{maxN} and \code{P_maxN} to find N corresponding to a distribution
#'		with a probability of \code{P_maxN} of exceeding \code{maxN}.}
#'	\item{\code{'lognormal'}}{Lognormal distribution with mean \code{p1}
#'		and std dev \code{p2}. Can alternatively specify \code{maxN} and 
#'		\code{P_maxN} to find the std dev corresponding to a distribution with 
#'		mean  = 0 and a probability of code{P_maxN} of exceeding \code{maxN}.}
#'	\item{\code{'poisson'}}{Poisson distribution with rate (lambda) = 
#'		\code{p1}. Can alternatively specify \code{maxN} and \code{P_maxN} to   
#'		find \code{p1} corresponding to a distribution with a probability of  
#'		\code{P_maxN} of exceeding \code{maxN}.}
#' Note that is \code{maxN} is specified, it will be used instead of 
#' \code{p1} or \code{p2}. In addition, continuous distributions are discretized
#' so that abundances can be used as measures of propagule supply.
#'
#' @note Users should take care when specifying correlations so that the 
#' results covariance matrix for niche parameters and global abundance is
#' positive definite. For example, if \code{r = 0.5} and global abundance is 
#' correlated with niche optima along environmental axis 1, it must also be 
#' correlated with niche optima along axis 2.
#'
#' @param N_S (required) number of species
#' @param nicheparms (required) a named list decribing the distributions from 
#' 	which niche optima and breadths are sampled. Contains:
#'	\decribe{
#'		\item{mu}{vector of length 2 with the maximum niche optima for each 
#'			environmental axis}
#'		\item{rho}{correlation between niche optima aross the two environmental
#'			axes}
#'		\item{sigma}{vector of length 2 with the means of the gamma 
#'			distributions from which niche breadths are sampled, corresponding
#'			to the two environmental axes}
#'		\item{alpha}{vector of length 2 with the shape parameters of the gamma
#'			distributions from which niche breadths are samples, corresponding
#'			to the two environmental axes}
#'		\item{r}{correlation between niche breadth across the two environmental
#'			axes}
#'	}
#' @param distribution (required) named list of parameters describing the 
#' form of the species abundance distribution. Contains:
#'	\decribe{
#'		\item{type}{character string indicating the form of the distribution.
#'			See details.}
#'		\item{maxN}{maximum abundance. Used when other distribution parameters
#'			are not supplied.}
#'		\item{P_maxN}{probability of sampling an abundance greater than 
#'			\code{maxN}. Used in conjunction with \code{maxN} to define the
#'			shape of the distribution when other parameters are not supplied.}
#'		\item{p1}{parameter controling the shape of the distribution. 
#'			Interpretation depends on \code{type}. See details.}
#'		\item{p2}{parameter controling the shape of the distribution. 
#'			Interpretation depends on \code{type}. See details.}
#'		\item{corr}{matrix containing correlations between global abundance and
#'			niche means and breadths. See details.

make_niches_gsad = function(N_S, nicheparms, distribution){

	# Create correlation matrix for generating mus, sigmas, and abundances
	S = diag(5)
	rownames(S) = c('mu1','mu2','sig1','sig2','abun')
	colnames(S) = c('mu1','mu2','sig1','sig2','abun')
	S['mu1','mu2'] = nicheparms$rho
	S['mu2','mu1'] = nicheparms$rho
	S['sig1','sig2'] = nicheparms$r
	S['sig2','sig1'] = nicheparms$r
	S[c('mu1','mu2','sig1','sig2'),'abun'] = distribution$corr
	S['abun', c('mu1','mu2','sig1','sig2')] = distribution$corr

	# Generate normal random variables with the given correlation structure
	x_norm = mvrnorm(n=N_S, mu=rep(0, 5), Sigma=S)

	# Transform to uniform on (0,1)
	U = pnorm(x_norm)

	# Niche optima
	# Transform to uniform on interval defined by nicheparms
	mu1 = qunif(U[,'mu1'], -nicheparms$mu[1], nicheparms$mu[1])
	mu2 = qunif(U[,'mu2'], -nicheparms$mu[2], nicheparms$mu[2])

	# Niche breadth
	# Transform to gamma with parameters defined by nicheparms
	theta = nicheparms$sigma/nicheparms$alpha # gamma distribution mean = scale*shape
	sig1 = qgamma(U[,'sig1'], shape=nicheparms$alpha[1], scale=theta[1])
	sig2 = qgamma(U[,'sig2'], shape=nicheparms$alpha[2], scale=theta[2])

	# Abundance
	if(distribution$type=='same') abuns = rep(1, N_S)
	
	if(distribution$type=='uniform'){
		if(is.null(distribution$maxN)){ 
			use_p = distribution$p1
		} else {
			use_p = distribution$maxN
		}
		abuns = qunif(U[,'abun'], 1, use_p)
	}
		
	if(distribution$type=='power'){
		if(!exists(distribution$maxN)){
			use_p = distribution$p1
		} else {
			
			# Stop if P_maxN not specified
			if(!exists(distribution$P_maxN)) stop('Must specify P_maxN.')
			
			use_p = 1/log(1/(distribution$maxN), base=distribution$P_maxN)
		}
		abuns = (1-U[,'abun'])^(-1/use_p)
	}
	if(distribution$type=='logseries'){
		# p2 = N, p1 = Fisher's alpha
		if(!exists(distribution$maxN)){
			use_N = distribution$p2
		} else {
			# Stop if P_maxN not specified
			if(!exists(distribution$P_maxN)) stop('Must specify P_maxN.')
		
			use_N = 1
			this_N = qls(1-distribution$P_maxN, use_N, distribution$p1)	
			while(this_N < distribution$maxN){
				use_N = use_N + 1
				this_N = qls(1-distribution$P_maxN, use_N, distribution$p1)	
			}
			if(this_N > 10) use_N = use_N - 1

		}
		abuns = qls(U[,'abun'], N = use_N, alpha = distribution$p1)		
	}
	if(distribution$type=='lognormal'){
		if(!exists(distribution$maxN)){
			use_mean = distribution$p1
			use_sd = distribution$p2
		} else {
			# Stop if P_maxN not specified
			if(!exists(distribution$P_maxN)) stop('Must specify P_maxN.')
			
			# Assume mean = 0
			use_mean = 0
			try_sig = 1
			this_N = qlnorm(1-distribution$P_maxN, use_mean, try_sig)
			i = 1
			while(abs(distribution$maxN - this_N) > 0.5){
				if(this_N > distribution$maxN){
					try_sig = try_sig - try_sig/2
					
				} else {
					try_sig = try_sig + try_sig/2
				}
				this_N = qlnorm(1-distribution$P_maxN, 0, try_sig)
			}
			use_sd = try_sig
		}

		# p1 = mean, p2 = sd
		abuns = ceiling(qlnorm(U[,'abun'], use_mean, use_sd))
	}
	if(distribution$type=='poisson'){
		# p1 = lambda
		if(!exists(distribution$maxN)){
			use_lamda = distribution$p1
		} else {
			# Stop if P_maxN not specified
			if(!exists(distribution$P_maxN)) stop('Must specify P_maxN.')
			
			use_maxN = distribution$maxN + 1
			try_lamda = use_maxN:1
			this_N = qpois(1-distribution$P_maxN, try_lamda)
			N_diffs = this_N - use_maxN
			i = 1
			while(!(0 %in% N_diffs)){
				new_start = which(abs(N_diffs)==min(abs(N_diffs)))
				try_lamda = seq(try_lamda[new_start+1], try_lamda[new_start-1], 10^(-i))
				this_N = qpois(1-distribution$P_maxN, try_lamda)
				N_diffs = this_N - use_maxN
				i = i + 1
			}
		
			use_lamda = min(try_lamda[which(N_diffs==0)])
		}

		abuns = qpois(U[,'abun'], use_lamda)+1
	}

	# Discretize so that it can be used in dbinom in the disperse function
	abuns = round(abuns, 0)
	
	# Return list of niches abuns gsad
	niches = array(c(mu1, sig1, mu2, sig2), dim=c(N_S, 2, 2), dimnames=list(1:N_S, c('mu','sigma'), 1:2))
	list(niches=niches, gsad=abuns)
}

