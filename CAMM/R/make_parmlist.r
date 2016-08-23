#' Create parameter list
#'
#' Creates a list of simulation parameters
#'
#' The function creastes a list of parameters required for running a 
#' simulation. Parameters are taken by default from the environment where
#' the function is called, or optionally from the environment provided
#' by \code{e}. Most parameters are \strong{required} for simulation and 
#' \emph{must be present in the environment or the function will fail}. 
#' Optional parameters are indicated.
#' 
#' Parameters for creating a metacommunity (see \code{\link{initialize_camm}}):
#' \itemize{
#'	\item \strong{\code{runID}} : name of this simulation run
#'	\item \strong{\code{S_a}} : number of host species
#'	\item \strong{\code{S_b}} : number of symbiont species
#'	\item \strong{\code{N_C}} : number of sites
#'	\item \strong{\code{N}} : number of microsites / individuals
#'	\item \strong{\code{rho_z}} : correlation between environmental axes
#'		across sites 
#'	\item \strong{\code{N_L}}: number of links in the association network 
#'		between host and symbiont species
#'	\item \strong{\code{topology}}: character string defining type of 
#'		association network topology
#'	\item \strong{\code{N_L}}: number of links in the association network 
#'		between host and symbiont species
#'	\item \strong{\code{comm_fill}}: (optional) whether initial metacommunity 
#'		should be \code{'empty'} or filled at \code{'random'}
#' }
#' Parameters for initializing a metacommunity with data 
#'(see \code{\link{initialize_camm}}, all parameters are optional):
#' \itemize{
#'	\item \strong{\code{topo_data}} : path to \code{csv} file with species
#'		association matrix
#'	\item \strong{\code{site_data}} : path to \code{csv} file with site
#'		environmental data
#'	\item \strong{\code{gsad_a_data}} : path to \code{txt} file with host
#'		global species abundances
#'	\item \strong{\code{gsad_b_data}} : path to \code{txt} file with symbiont
#'		global species abundances
#' }
#' Parameters for creating species niches (see \code{\link{make_niches_gsad}}):
#' \itemize{
#'	\item \strong{\code{mu_a1}} and \strong{\code{mu_a2}} : maximum niche
#'		optima for host species along environmental axes 1 and 2
#'	\item \strong{\code{mu_b1}} and \strong{\code{mu_b2}} : maximum niche
#'		optima for symbiont species along environmental axes 1 and 2
#'	\item \strong{\code{rho_a}} and \strong{\code{rho_b}}: correlation of niche
#'		optima between environmental axes across host and symbiont species
#'	\item \strong{\code{sigma_a1}} and \strong{\code{sigma_a2}} : expected
#'		niche breadth of host species along environmental axes 1 and 2
#'	\item \strong{\code{sigma_b1}} and \strong{\code{sigma_b2}} : expected
#'		niche breadth of symbiont species along environmental axes 1 and 2
#'	\item \strong{\code{alpha_a1}} and \strong{\code{alpha_a2}} : shape
#'		parameter of the gamma distribution from which host niche breadths
#'		are drawn for environmental axes 1 and 2
#'	\item \strong{\code{alpha_b1}} and \strong{\code{alpha_b2}} : shape
#'		parameter of the gamma distribution from which host niche breadths
#'		are drawn for environmental axes 1 and 2
#'	\item \strong{\code{r_a}} and \strong{\code{r_b}}: correlation of niche
#'		breadths between environmental axes across host and symbiont species
#'	\item \strong{\code{gsad_dist_a}} and \strong{\code{gsad_dist_b}}: list
#'		of paramters describing the global species abundance distribution for
#'		hosts ans symbionts
#' }
#' Parameters for running the simulation:
#' \itemize{
#'	\item \strong{\code{omega}}: degree of mutualism: \code{0} = no 
#'		facilitation, \code{1} = obligate mutualism 
#'		(see \code{\link{calc_probs}})
#'	\item \strong{\code{assoc_prob_func}}: function that defines the
#'		probability that a host and symbiont will form an association, given
#'		that they have the capacity to. Should accept three arguments: 
#'		\code{env} (vector of environmental conditions), \code{niche_a} 
#'		(matrix defining the host's niche) and \code{niche_b} (matrix defining
#'		the symbiont's niche).
#'	\item \strong{\code{mort_rate}}: probability that an established association
#'		dies in a timestep
#'	\item \strong{\code{mort_rate_a}} and \strong{\code{mort_rate_b}}: local
#'		extinction rate of unassociated hosts and symbionts from a species pool
#'		relative to the mortality rate of established associations (e.g., a
#'		value of \code{10} will cause cause unassociated mutualists to be removed
#'		from a species pool at \code{10x} the rate of associated mutualists)
#' }
#'
#' @param e environment where parameters are located
#' @return a named list of paramters
#'
#' @export
make_parmlist = function(e=parent.frame()){
	parms = list(
		runID = e$runID,
#		a_name = e$a_name, 
#		b_name = e$b_name,
		S_a = e$S_a,
		S_b = e$S_b,
		N_C = e$N_C,
		N = e$N,
		rho_z = e$rho_z, 
		mu_a1 = e$mu_a1,
		mu_a2 = e$mu_a2,
		rho_a = e$rho_a,
		mu_b1 = e$mu_b1,
		mu_b2 = e$mu_b2,
		rho_b = e$rho_b, 
		sigma_a1 = e$sigma_a1,
		sigma_a2 = e$sigma_a2,
		alpha_a1 = e$alpha_a1,
		alpha_a2 = e$alpha_a2, 
		r_a = e$r_a,
		sigma_b1 = e$sigma_b1, 
		sigma_b2 = e$sigma_b2,
		alpha_b1 = e$alpha_b1,
		alpha_b2 = e$alpha_b2,
		r_b = e$r_b,	
		gsad_dist_a = e$gsad_dist_a,
		gsad_dist_b = e$gsad_dist_b,
		N_L = e$N_L,
		topology = e$topology,
		omega = e$omega,
		assoc_prob_func = e$assoc_prob_func,
		mort_rate = e$mort_rate,
		mort_rate_a = e$mort_rate_a,
		mort_rate_b = e$mort_rate_b
	)
	
	# Add optional parameters
	if(exists('comm_fill', e)) parms = c(parms, comm_fill = e$comm_fill)
	if(exists('topo_data', e)) parms = c(parms, topo_data = e$topo_data)
	if(exists('site_data', e)) parms = c(parms, site_data = e$site_data)
	if(exists('gsad_a_data', e)) parms = c(parms, gsad_a_data = e$gsad_a_data)
	if(exists('gsad_b_data', e)) parms = c(parms, gsad_b_data = e$gsad_b_data)
	
	# Return list
	parms
	
}
