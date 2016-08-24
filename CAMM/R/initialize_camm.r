### TO DO: MODIFY FUNCTION TO ACCEPT EXISTING METACOMM OBJECTS
#' Initialize simulation
#'
#' Initializes a simulation for a given set of parameters
#'
#' The function creates host and symbiont species pools, global species 
#' abundance distributions (gsad), an association network (and named links),
#' an association probability matrix, a set of sites with environmental 
#' conditions, and an empty metacommunity at those sites. Parameters for
#' creating these objects are read from a file (\code{parm_file}) or from 
#' the environment where the function is called. Required parameters are:
#' \describe{
#'	\item{\code{S_a} and \code{S_b}}{number of host and symbiont species}
#'	\item{\code{N_L}}{number of links in the association network} 
#'	\item{\code{topology}}{association network topology (\code{'one2one'},
#'		\code{'one2many'} or \code{'many2many'}). See \code{\link{make_topo}}.}
#'	\item{\code{N_C}}{number of sites}
#'	\item{\code{N}}{number of microsites/individuals per site}
#'	\item{\code{rho_z}}{correlation between environmental axes across sites}
#'	\item{\code{niche_parms_a} and \code{niche_parms_b}}{named list of parameters
#'		used to create host (a) and symbiont (b) niches. 
#'		See \code{\link{make_niches_gsad}}.}
#'	\item{\code{gsad_dist_a} and \code{gsad_dist_b}}{named list of parameters
#'		used to create host (a) and symbiont (b) global species abundance 
#'		distributions. See \code{\link{make_niches_gsad}}.}
#' }
#' Optional parameters are:
#' \describe{
#'	\item{comm_fill}{character string indicating how initial metacommunity
#'		should be filled. Defaults to \code{'empty'}. See \code{link{make_comm}}.}
#'	\item{\code{topo_data}}{path to \code{csv} file with species
#'		association matrix. Path should be relative to the location of 
#'		\code{parm_file}, if provided. File should contain a binary matrix  
#'		without row or column names. If provided, overrides values for  
#'		\code{S_a}, \code{S_b}, \code{N_L}, and \code{topology}.}
#'	\item{\code{site_data}}{path to \code{csv} file with site
#'		environmental data. Path should be relative to the location of 
#'		\code{parm_file}, if provided. File should contain two columns: one 
#'		for each environmental axis. If provided, overrides values for 
#'		\code{N_C} and \code{rho_z}.}
#'	\item{\code{gsad_a_data}}{path to \code{txt} file with host
#'		global species abundances. Path should be relative to the location of 
#'		\code{parm_file}, if provided. File should contain one line with species
#'		abundances (as integers) separated by white space. If provided, overrides 
#'		values for \code{gsad_dist_a} and will not be correlated with niches.}
#'	\item{\code{gsad_b_data}}{path to \code{txt} file with symbiont
#'		global species abundances. Path should be relative to the location of 
#'		\code{parm_file}, if provided. File should contain one line with species
#'		abundances (as integers) separated by white space. If provided, overrides 
#'		values for \code{gsad_dist_b} and will not be correlated with niches.}
#' }
#'
#' @param parm_file path to file with parameters
#' @param save_start logical indicating whether metacommunity objects should be
#' 	saved to an .RData file. Defaults to \code{FALSE}.
#' @param runID character string identifying this simulation. Required if 
#' 	\code{save_start=TRUE}.
#' @param save_dir directory where simulation objects should be saved. Defaults
#' 	to current directory.
#' @return a named list containing all objects needed to run a simulation
#'
#' @seealso \code{\link{make_topo}}, \code{\link{make_niches_gsad}}, 
#' 	\code{\link{make_comm}}, \code{\link{make_sites}}
#'
#' @export

initialize_camm = function(parm_file=NA, save_start = F, runID=NA, save_dir='./'){

	# Read in parameter values
	if(!is.na(parm_file)) source(parm_file)
	
	# Get directory name of parameter file
	parm_dir = ifelse(!is.na(parm_file), dirname(parm_file), getwd())

	# Assign optional parameters if not present
	if(!exists('comm_fill', parent.frame())) comm_fill='empty'
	
	## Instantiate communities and mutualistic network
	
	# If data for association network is provided, read it in, otherwise create a random network according to parameters
	if(exists('topo_data', parent.frame())){
		# location of TOPO data is relative to the location of the parameter file
		oldwd = getwd()
		setwd(parm_dir)
		topo = tryCatch(read.csv(topo_data, header=F), error = function(e){
			setwd(oldwd)
			stop(paste('Could not load network data from', topo_data))
		})
		setwd(oldwd)
		S_a = nrow(topo)
		S_b = ncol(topo)
		N_L = sum(topo>0)
		topoA = ifelse(sum(rowSums(topo)>1)>0, 'many', 'one')
		topoB = ifelse(sum(colSums(topo)>1)>0, 'many', 'one')
		topology = paste0(topoA, '2', topoB)
	} else {
		topo = make_topo(S_a, S_b, N_L, topology) # S_a x S_b matrix of association probabilities
	}
	
	# Create names for species associations
	topo_names = name_topo(topo) # S_a x S_b matrix with integers labeling specific associations
	
	# If data for sites is provided, read it in, otherwise create random sites
	if(exists('site_data', parent.frame())){
		# location of site data is relative to the location of the parameter file
		oldwd = getwd()
		setwd(parm_dir)
		sites = tryCatch(read.csv(site_data), error = function(e){
				setwd(oldwd)
				stop(paste('Could not load site data from', site_data))
		})
		setwd(oldwd)
		N_C = nrow(sites)
		rho_z = cor(sites[,1], sites[,2])
	} else {
		sites = make_sites(N_C, rho_z) # N_C x 2 matrix of env values
	}

	
	# Fill initial community
	comm = make_comm(N_C, N, N_L, fill=comm_fill) # initial N_C x N matrix of integers indicating which association is present
	poolA = calc_pool(comm, topo_names, 1) # N_C x S_a matrix of presence of mutualist a in communities based on associations present
	poolB = calc_pool(comm, topo_names, 2) # N_C x S_b matrix of presence of mutualist b in communities based on associations present

	
	# Generate random niches for each mutualist
	niche_gsad_a = make_niches_gsad(S_a, nicheparms_a, gsad_dist_a)
	niche_gsad_b = make_niches_gsad(S_b, nicheparms_b, gsad_dist_b)
	niches_a = niche_gsad_a$niches # array of S_a x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
	niches_b = niche_gsad_b$niches # array of S_b x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)

	
	# If data on gsads provided use this, otherwise use randomly generated gsads
	if(exists('gsad_a_data', parent.frame())){
		# location of gsad data is relative to the location of the parameter file
		oldwd = getwd()
		setwd(parm_dir)
		gsad_a = tryCatch(as.integer(read.table(gsad_a_data)), error = function(e){
			setwd(oldwd)
			stop(paste('Could not read abundances from', gsad_a_data))
		})
		setwd(oldwd)
		if(length(gsad_a)!=S_a) stop(paste('Supplied GSAD does not have', S_a, 'species.'))
	} else {
		gsad_a = niche_gsad_a$gsad
	}
	if(exists('gsad_b_data', parent.frame())){
		# location of gsad data is relative to the location of the parameter file
		oldwd = getwd()
		setwd(parm_dir)
		gsad_b = tryCatch(as.integer(read.table(gsad_b_data)), error = function(e){
			setwd(oldwd)
			stop(paste('Could not read abundances from', gsad_b_data))
		})
		setwd(oldwd)
		if(length(gsad_a)!=S_b) stop(paste('Supplied GSAD does not have', S_b, 'species.'))
	} else {
		gsad_b = niche_gsad_b$gsad
	}


	# Calculate N_C x N_L matrix of association probabilities at each site
	# These govern the probability that partner A will establish 
	#  given that its mutualist B is present and weighted by the strength of mutualism (omega)
	assoc_probs = matrix(0, N_C, N_L)
	for(i in 1:N_C){
	for(j in 1:N_L){
		env = sites[i,]

		partners = which(topo_names==j, arr.ind=T)
		n_a = niches_a[partners[1],,]
		n_b = niches_b[partners[2],,]

		assoc_probs[i,j] = assoc_prob_func(env, n_a, n_b)
	}}

	# Save community parameters
	if(save_start) save(topo, topo_names, sites, niches_a, niches_b, gsad_a, gsad_b, assoc_probs, file=file.path(save_dir, paste0(runID, '_simobject.RData')))
	
	# Return initial metacommunity
	metacomm = list(topo=topo, topo_names=topo_names, niches_a=niches_a, niches_b=niches_b, 
		gsad_a=gsad_a, gsad_b=gsad_b, assoc_probs=assoc_probs, sites=sites, comm=comm, poolA=poolA, poolB=poolB)
	metacomm
}