## This script holds functions used by the community assembly of mutualists model (CAMM)
library(vegan) # vegdist
require(poweRlaw)
require(sads)
require(permute)
require(MASS) # mvrnorm

### Functions for Initializing Simulation ###

# A function that generates a random bipartite graph between sets A and B
# under the following constraints:
# 	S_a = size of A
# 	S_b = size of B
# 	N_L = number of edges
#	topology = 'one2one', 'one2many', 'many2many'
# Returns a S_a x S_b binary matrix

make_topo = function(S_a, S_b, N_L, topology){
	# Catch errors
	if(N_L < S_a | N_L < S_b) stop('Cannot find partners for all mutualists when N_L < S_a or S_b')
	if(topology=='one2one' & S_a!=S_b) stop('Cannot form one-to-one topology when S_a != S_b')
	if(topology=='one2one' & S_a!=N_L) stop('Cannot form one-to-one topology when S_a != N_L')
	if(topology=='one2many' & N_L <= S_b) stop('Cannot form one-to-many topology when N_L <= S_b')
	if(topology!='many2many' & N_L > S_a) stop('Topology must be many-to-many if N_L > S_a')
	if(N_L > S_a*S_b) stop(paste('Cannot form', N_L, 'associations with', S_a,'and', S_b, 'partners'))

	# Make empty association matrix
	topo = matrix(0, nrow=S_a, ncol=S_b)	

	# Choose one partner for each mutualist 
	partners_a = sample(S_a, S_a, replace=F)
	partners_b = sample(S_b, S_b, replace=F)

	# If topology is many2many add extra links for mutualist a
	if(topology=='many2many') partners_a = c(partners_a, sample(S_a, N_L-S_a, replace=T))

	# If topology is one2many or many2many add extra linkes for mutualist b
	if(topology %in% c('one2many','many2many')) partners_b = c(partners_b, sample(S_b, N_L-S_b, replace=T))

	# Shuffle order
	partners_a = sample(partners_a, replace=F)
	partners_b = sample(partners_b, replace=F)
	

	# Add links to network
	for(i in 1:N_L){
		topo[partners_a[i],partners_b[i]] = 1
	}

	# If there are duplicated links add more random links until at N_L.
	add_links = N_L - sum(topo)
	while(add_links > 0){
		add_a = sample(S_a, add_links, replace=T)
		add_b = sample(S_b, add_links, replace=T)
		for(i in 1:add_links){
			topo[add_a[i], add_b[i]] = 1
		}
		add_links = N_L - sum(topo)
	}

	topo
}

# A function that labels each edge of a bipartite graph with a unique integer
#	topo = a binary matrix with 1 indicating association
# Returns a matrix with the same structure but with 1 replaced by integer labels

name_topo = function(topo){
	topo_labeled = t(topo)
	N_L = sum(topo_labeled==1)
	topo_labeled[topo_labeled==1] = 1:N_L

	t(topo_labeled)
}


# A function that generate 2 random environmental variables for:
# 	N_C = number of sites
#	rho_z = correlation between variables
# Returns an N_C x 2 matrix of values drawn from a multivariate gaussian distribution with mean 0 and std deviation 1.

make_sites = function(N_C, rho_z){
	sd = c(1,1)
	mu = c(0,0)

	# Define covariance matrix
	S = matrix(c(sd[1],rho_z,rho_z,sd[2]), nrow=2, ncol=2)

	# Generate correlated gaussian rv using MASS package
	sites = mvrnorm(n=N_C, mu=mu, Sigma=S)
	
	# Return sites
	sites
}


# A function that generates an initial set of communities
# Currently generates empty communities
#	N_C = number of communities
#	N = number of individuals in each community
#	N_S = number of potential species (associations)

make_comm = function(N_C, N, N_S){
	comm = matrix(0, nrow=N_C, ncol=N)
	comm
}

# A function that generates a random community for the purpose of testing the simulation
#	N_C = number of communities
#	N = number of individuals in each community
#	N_S = number of potential species

rand_comm = function(N_C, N, N_S){
	comm = matrix(sample(0:N_S, N_C*N, replace=T), N_C, N)	
}



# A function that generates random 2D gaussian niches.
# Niche optima (mu) are drawn from a uniform distribution and niche width (sigma) from a gamma distribution
#	N_S = number of species
#	nicheparms = a list of parameters controlling the shape of the distribution from which niches are sampled
#		mu = a vector of length 2 with the maximum niche optima
#		rho = the correlation between the niche optima
#		sigma = a vector of length 2 with the means of the gamma distributions
#		alpha = a vector of length 2 with the shape parameter of the gamma distributions
#		r = the correlation between niche widths

make_niches = function(N_S, nicheparms){
	
	# Generate niche optima as 2 correlated gaussian variables using MASS package
	S = matrix(c(1,nicheparms$rho,nicheparms$rho,1),2,2)
	mus_norm = mvrnorm(n=N_S, mu=c(0,0), Sigma=S)
	
	# Transform to uniform on (0,1)
	U = pnorm(mus_norm)

	# Transform to uniform on interval defined by nicheparms
	mu1 = qunif(U[,1], -nicheparms$mu[1], nicheparms$mu[1])
	mu2 = qunif(U[,2], -nicheparms$mu[2], nicheparms$mu[2])

	# Generate niche breadths as 2 correlated gaussian variables
	S = matrix(c(1, nicheparms$r, nicheparms$r, 1),2,2)
	sigmas_norm = mvrnorm(n=N_S, mu=c(0,0), Sigma=S)

	# Transform to uniform on (0,1)
	U = pnorm(sigmas_norm)	
	
	# Transform to gamma with parameters defined by nicheparms
	theta = nicheparms$sigma/nicheparms$alpha # gamma distribution mean = scale*shape
	sig1 = qgamma(U[,1], shape=nicheparms$alpha[1], scale=theta[1])
	sig2 = qgamma(U[,2], shape=nicheparms$alpha[2], scale=theta[2])

	# Return array of niches
	niches = array(c(mu1, sig1, mu2, sig2), dim=c(N_S, 2, 2), dimnames=list(1:N_S, c('mu','sigma'), 1:2))
	niches
}

# THIS FUNCTION SHOULD NOT BE USED ANY LONGER
# A function that generates a global species abundance distribution.
# Species abundances may be correlated with niche optima or niche breadth or mutualist breadth.
# S = number of species
# distribution = list of parameters describing the functional form of the SAD
make_gsad = function(S, distribution){
	
	if(distribution$type=='same'){
		abuns = rep(1, S)
	}
	if(distribution$type=='uniform'){
		abuns = runif(S)
	}
	if(distribution$type=='power'){
		# p1 = alpha
		abuns = rpldis(S, 1, distribution$p1)
	}
	if(distribution$type=='logseries'){
		N = 1000000
		
		# p1 = N/alpha
		alpha = N / distribution$p1

		# Find the n below which 99.9% of species abundances occur.
		n = qls(0.999, N, alpha)
		
		abuns = sample(1:n, S, prob=dls(1:n, N, alpha), replace=T)		
	}
	if(distribution$type=='lognormal'){	
		# p1 = mean, p2 = sd
		abuns = rlnorm(S, distribution$p1, distribution$p2)
	}
	if(distribution$type=='poisson'){
		# p1 = lambda
		abuns = rpois(S, distribution$p1)+1
	}
	
	abuns	
}

# A function that concurrently generates random 2D gaussian niches and a species abundance distribution that is correlated with niche parameters.
# Niche optima (mu) are drawn from a uniform distribution and niche width (sigma) from a gamma distribution
# NOTE: be careful with choice of covariances so that covariance matrix works (positive definite)
#	N_S = number of species
#	nicheparms = a list of parameters controlling the shape of the distribution from which niches are sampled
#		mu = a vector of length 2 with the maximum niche optima
#		rho = the correlation between the niche optima
#		sigma = a vector of length 2 with the means of the gamma distributions
#		alpha = a vector of length 2 with the shape parameter of the gamma distributions
#		r = the correlation between niche widths
# 	distribution = list of parameters describing the functional form of the SAD
#		type = vector indicating the form of the SAD
#		maxN = maximum abundance, used to find SAD parameters if not supplied
#		P_maxN = probability of getting an abundance greater than maxN
#		p1 = parameter controling the shape of the SAD (only certain types)
#		p2 = parameter controling the shape of the SAD (only certain types)
#		corr = matrix giving correlations between global abundance and niche means and breadths (rows are env axes)
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
		if(is.null(distribution$maxN)){
			use_p = distribution$p1
		} else {
			use_p = 1/log(1/(distribution$maxN), base=distribution$P_maxN)
		}
		abuns = (1-U[,'abun'])^(-1/use_p)
	}
	if(distribution$type=='logseries'){
		# p2 = N, p1 = Fisher's alpha
		if(is.null(distribution$maxN)){
			use_N = distribution$p2
		} else {
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
		if(is.null(distribution$maxN)){
			use_mean = distribution$p1
			use_sd = distribution$p2
		} else {
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
		if(is.null(distribution$maxN)){
			use_lamda = distribution$p1
		} else {
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



# A function that initializes a simulation run with a given set of parameters
# If parm_file is missing, then paramters are taken from the environment where the function is called
# parm_file = text file with parameter values
# save_start = flag indicating whether niches and topology should be saved
# runID = text string indicating how initial settings should be saved
initialize_camm = function(parm_file=NA, save_start = F, runID=NA, save_dir='./'){
	# Read in parameter values
	if(!is.na(parm_file)) source(parm_file)
	
	# Instantiate communities and mutualistic network
	topo = make_topo(S_a, S_b, N_L, topology) # S_a x S_b matrix of association probabilities
	topo_names = name_topo(topo) # S_a x S_b matrix with integers labeling specific associations
	sites = make_sites(N_C, rho_z) # N_C x 2 matrix of env values
	comm = make_comm(N_C, N, N_L) # initial N_C x N matrix of integers indicating which association is present
	poolA = calc_pool(comm, topo_names, 1) # N_C x S_a matrix of presence of mutualist a in communities based on associations present
	poolB = calc_pool(comm, topo_names, 2) # N_C x S_b matrix of presence of mutualist b in communities based on associations present

	# Generate random niches for each mutualist
	niche_gsad_a = make_niches_gsad(S_a, nicheparms_a, gsad_dist_a)
	niche_gsad_b = make_niches_gsad(S_b, nicheparms_b, gsad_dist_b)
	niches_a = niche_gsad_a$niches # array of S_a x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
	niches_b = niche_gsad_b$niches # array of S_b x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
	gsad_a = niche_gsad_a$gsad
	gsad_b = niche_gsad_b$gsad

	# Calculate N_C x N_L matrix of association probabilities at each site
	# These govern the probability that partner A will establish 
	#  given that its mutualist B is present and weighted by the strength of mutualism (omega)
	assoc_probs = matrix(0, N_C, N_L)
	for(i in 1:N_C){
	for(j in 1:N_L){
		env = sites[i,]

		partners = which(topo_names==j, arr.ind=T)
		n_a = niches_a[partners[1],'mu',]
		n_b = niches_b[partners[2],'mu',]

		assoc_probs[i,j] = assoc_prob_func(env, n_a, n_b)
	}}

	# Save community parameters
	if(save_start) save(topo, topo_names, sites, niches_a, niches_b, gsad_a, gsad_b, assoc_probs, file=paste0(save_dir, runID, '_simobject.RData'))
	
	# Return initial metacommunity
	metacomm = list(topo=topo, topo_names=topo_names, niches_a=niches_a, niches_b=niches_b, 
		gsad_a=gsad_a, gsad_b=gsad_b, assoc_probs=assoc_probs, sites=sites, comm=comm, poolA=poolA, poolB=poolB)
	metacomm
}

### Functions for Running Simulation ###

# A function that calculates the probability of survival in environment x given niche parameters
# If x is a vector, probabilities are multiplied.
#	x = vector of environmental values
#	mu = vector of niche optima
#	sigma = vector of niche breadths
niche_func = function(x, mu, sigma){
	N = length(x)
	pvals = sapply(1:N, function(i) 1/exp((x[i]-mu[i])^2/(2*sigma[i]^2)))

	prod(pvals)
}

# A function that surveys which mutualist partners are present in a community
#	comm = matrix indicating which association is present in each position of a community
#	topo_names = bipartite mutualist network with integers labeling associations
#	partner = integer indicating whether this is for mutualist 'a' (1) or 'b' (2)

calc_pool = function(comm, topo_names, partner){
	# Find number of communities and species
	N_C = nrow(comm)
	N_S = dim(topo_names)[partner]

	# For each row of comm (representing a community) find the associations present and find the partners present
	pool = matrix(0, nrow=N_C, ncol=N_S)
	for(i in 1:N_C){

		# Find the associations present
		assocs = unique(comm[i,])
		assocs = assocs[assocs!=0]

		# Find the partners present
		for(x in assocs){
			j = which(topo_names==x, arr.ind=T)[,partner]

			# Assign 1 for presence in community i
			pool[i,j] = 1
		}
	}

	# Return partner pool
	pool
}


# A function that causes random mortality of mutualists in a community. Mutualists involved in an association cannot die.
#	comm = matrix indicating which association is present in each position of a community
#	topo_names = bipartite mutualist network with integers labeling associations
#	pool = site X mutualist binary matrix indicating which mutualist species are present
# 	mortality = probability that an unassociated mutualist dies
#	partner = integer indicated which partner this applies to
die = function(comm, topo_names, pool, mortality, partner){
	
	# Find pool of mutualists that are not in an association and can die
	comm_pool = calc_pool(comm, topo_names, partner)
	not_associated = which(pool!=comm_pool)

	# For each unassociated partner in each community do random mortality
	survive = runif(length(not_associated)) > mortality
	pool[not_associated] = pool[not_associated]*survive

	# Return new pool
	pool
}

# A function that causes random dispersal of mutualists into a community.
# For now, colonization events independent with probability proportional to the global species abundance distribution for each partner.
#	sites = a 2 column matrix with the environmental conditions at each site
#	niches = a 2 column matrix giving niche optima and breadths for each mutualist species
#	pool = site X mutualist binary matrix indicating which mutualist species are present

disperse = function(sites, niches, pool, gsad){
	
	# Find number of communities and species
	N_C = nrow(pool)
	N_S = ncol(pool)

	# For each species and each site calculate the probability of establishment based on joint gaussian niches on each environmental variable
	probs = matrix(0,N_C, N_S)	

	for(i in 1:N_C){
	for(j in 1:N_S){
		p = niche_func(sites[i,], niches[j,'mu',], niches[j,'sigma',])
		probs[i,j] = 1 - dbinom(0, gsad[j], p) # Probability of at least one propagule establishing with propagule pressure equal to global abundance
	}}

	# Generate random dispersal
	rands = runif(N_C*N_S) 
	emmigrants = rands <= probs
	
	# Add emmigrants to existing pool
	pool = pool + emmigrants > 0
	pool = matrix(as.numeric(pool), N_C, N_S)

	# Return new pool of partners
	pool	
}


# A function that calculates the transition probability matrix for each site
# Currently species cannot displace one another
#	sites = a 2 column matrix with the environmental conditions at each site
#	niches_a = a 2 column matrix giving niche optima and breadths for each mutualist species a
#	niches_b = a 2 column matrix giving niche optima and breadths for each mutualist species b
#	topo_names = topo_names = bipartite mutualist network with integers labeling associations
#	poolA = site X mutualist binary matrix indicating which mutualist species 'a' are present
#	poolB = site X mutualist binary matrix indicating which mutualist species 'b' are present
#	mort_rate =  probability that an association dies
#	assoc_prob = an N_C x N_L matrix giving the probability that an association forms in a given environment if both partners are present

calc_probs = function(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs, omega){
	# Calculate number of associations and communities
	S_a = nrow(topo_names)
	S_b = ncol(topo_names)
	N_L = max(topo_names)
	N_C = nrow(sites)

	# Create empty array that holds transition matrices for each site
	Tarr = array(NA, dim=c(N_L+1, N_L+1, N_C), dimnames=list(0:N_L, 0:N_L, 1:N_C))

	# A 2 x N_C x N_L matrix indicating whether partner A or B is present in each community
	AB_present = sapply(1:N_L, function(x){
		partners = which(topo_names==x, arr.ind=T)
		rbind(poolA[,partners[1]],poolB[,partners[2]])
	}, simplify='array')	


	# Calculate probabilities of establishment 0 -> i
	establish = apply(AB_present, c(2,3), function(x){
		# If partner A is absent P = 0
		# Otherwise, if partner B is absent, P = 1-omega
		# Or, if partner B is also present, P = 1
		ifelse(x[1]==0, 0, ifelse(x[2]==0, 1-omega, 1))
	})*assoc_probs 
	
	# Calculate probability that no species establishes
	no_establish = apply(establish, 1, function(p) prod(1-p))

	# Scale each P(0 -> i) so that rows sum to 1. A species has P=establish only if it is the only species that can establish. 
	establish = establish/ifelse(rowSums(establish)>0, rowSums(establish),1)*(1-no_establish)# scale establishment probabilities so that 0 
	Tarr[1,,] = t(cbind(no_establish, establish))
	
	# Add probability of death
	# For now, constant for any established species
	Tarr[2:(N_L+1),1,] = mort_rate	
	
	# For each species add transition probabilities
	for(i in 1:N_L){
		# Probability of competitive displacement
		# For now, no competitive displacement
		# With competition, entries should be the product of:
			# probability that the competitor can establish 
			# a (environmentally dependent) competitive  ability coefficient
			# scaled to sum=1 if sum > 1
		compete = matrix(0, N_C, N_L-1)
		colnames(compete) = (1:N_L)[1:N_L != i]
		compete = (1-Tarr[i+1,1,])*compete # Multiply competitive probabilities by P(i survives)

		# Probability of remaining
		remain = 1 - Tarr[i+1,1,]- rowSums(compete)

		Tarr[i+1,i+1,] = remain
		Tarr[i+1, colnames(compete),] = t(compete)		

	}
	# Check that rows sum to 1
	if(sum(round(apply(Tarr, c(1,3), sum),10)!=1)>0) stop('Transition matrix rows do not sum to 1')

	# Return transition array
	Tarr
}


# A function that randomly transitions a place in a community based on transition probabilities
#	probs = transition probability vector

transition = function(probs){
	as.numeric(sample(names(probs), 1, prob=probs))
}

# A function for running a simulation that has already been initialized
# Can be run in the environment that contains the simulation objects or these objects can be passed in a list called metacomm
# Note: other parameter (mort_rate, omega) are taken from the environment where the function is called
# Run simulation for set number of time steps ('fixed') or until set number of chains converge ('converge')
# metacomm = named list describing a metacommunity (such as returned by initialize_camm function)
#	must contain: sites, comm, poolA, poolB, niches_a, niches_b, topo, topo_names, assoc_probs, gsad_a, gsad_b
# sim_mode = 'converge' or 'fixed'
# reps = If mode is 'fixed', must supply number of time steps
# nchains = If mode is 'converge', must supply number of chains (replicate simulation runs used to check convergence)
# save_steps = If absent, all steps are saved in an array. Otherwise, a vector indicating which steps to save in an array. Simulations run in mode 'converge' must have save_steps=NA
run_camm = function(metacomm=NULL, sim_mode, reps=NA, nchains=NA, save_steps=NA){
	if(length(metacomm)>0){
		comm = metacomm$comm
		poolA = metacomm$poolA
		poolB = metacomm$poolB
		sites = metacomm$sites
		assoc_probs = metacomm$assoc_probs
		gsad_a = metacomm$gsad_a
		gsad_b = metacomm$gsad_b
		niches_a = metacomm$niches_a
		niches_b = metacomm$niches_b
		topo = metacomm$topo
		topo_names = metacomm$topo_names	
	}

	if(sim_mode=='fixed'){
	
	# Make vector of steps to save that includes timepoints in reps (if it is a vector)
	if(is.na(save_steps)) save_steps = 0:max(reps)
	save_steps = unique(c(save_steps, reps))
	save_steps = save_steps[order(save_steps)]

	# Generate empty array to hold changes in community and add initial community
	nsteps = length(save_steps)
	
	comm_records = array(NA, c(N_C, N, nsteps))
	poolA_records = array(NA, c(N_C, S_a, nsteps)) 
	poolB_records = array(NA, c(N_C, S_b, nsteps))
	Tmat_records = array(NA, c(N_C, N_L+1, N_L+1, nsteps))

	if(is.na(save_steps)|(0 %in% save_steps)){
		comm_records[,,1] = comm
		poolA_records[,,1] = poolA
		poolB_records[,,1] = poolB
	}

	# Run simulation
	for(step in 1:reps){
		# Mutualist mortality: mutualists present as associations cannot die
		poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
		poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)

		# Mutualists disperse into communities independently and with probability based on niche-based filtering	
		poolA = disperse(sites, niches_a, poolA, gsad_a)
		poolB = disperse(sites, niches_b, poolB, gsad_b)
	
		# Calculate transition matrices for each site
		T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs, omega) # array of square matrices of transition probabilities from one association to another

		# Identify current state of each space and transition based on random numbers
		new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
		for(i in 1:nrow(comm)){
		for(j in 1:ncol(comm)){
			new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
		}}
		comm = new_comm

		# Save communities
		if(step %in% save_steps){
			i = which(save_steps==step)
			comm_records[,,i] = comm
			poolA_records[,,i] = poolA
			poolB_records[,,i] = poolB
			Tmat_records[,,,i] = T_mat
		}
	}
	} # CLOSES if(sim_mode=='fixed') 

	# Run simulation until N chains converge on equilibrial community dynamics
	# Currently uses Rhat statistic from Brooks and Gelman 1998
	# This compares the within chain percentile interval to the interval calculated across all chains
	if(sim_mode=='converge'){

	# Split community into nchain identical instances
	comm_arr = array(comm, dim=c(nrow(comm), ncol(comm), nchains))
	poolA_arr = array(poolA, dim=c(nrow(poolA), ncol(poolA), nchains))
	poolB_arr = array(poolB, dim=c(nrow(poolB), ncol(poolB), nchains))
	Tmat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs)
	Tmat_arr = array(Tmat, dim=c(dim(Tmat), nchains))
	Rhat_arr = array()

	# Generate empty array to hold changes in community and add initial community
	comm_records = array(comm_arr, dim=c(dim(comm_arr), 1)) 
	poolA_records = array(poolA_arr, dim=c(dim(poolA_arr),1)) 
	poolB_records = array(poolB_arr, dim=c(dim(poolB_arr),1))
	Tmat_records = array(Tmat_arr, dim=c(dim(Tmat_arr),1))

	# Run simulation until convergence
	window_size = 50 # Number of iterations to run before re-calculating convergence criteria
	burnin = 100
	thin = 1
	Rhat_tol = 0.1 # Tolerance for R_hat statistic to be considered converged

	# Set initial counters
	converged = F
	step = 1

	# COMMUNITIES ARE CONVERGING VERY SLOWLY UNDER THESE CRITERIA- MAY WANT TO RECONSIDER
	while(!converged){
		# Do one step through the simulation
		for(k in 1:nchains){
			comm = comm_arr[,,k]
			poolA = poolA_arr[,,k]
			poolB = poolB_arr[,,k]

			# Mutualist mortality: mutualists present as associations cannot die
			poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
			poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)

			# Mutualists disperse into communities independently and with probability based on niche-based filtering	
			poolA = disperse(sites, niches_a, poolA, gsad_a)
			poolB = disperse(sites, niches_b, poolB, gsad_b)

			# Calculate transition matrices for each site
			T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs, omega) # array of square matrices of transition probabilities from one association to another
			Tmat_arr[,,,k] = T_mat

			# Identify current state of each space and transition based on random numbers
			new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
			for(i in 1:nrow(comm)){
			for(j in 1:ncol(comm)){
				new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
			}}
			comm = new_comm

			# Save communities
			comm_arr[,,k] = comm
			poolA_arr[,,k] = poolA
			poolB_arr[,,k] = poolB
		}
	
		# Save communities from this step
		Tmat_records = abind(Tmat_records, Tmat_arr)
		comm_records = abind(comm_records, comm_arr)
		poolA_records = abind(poolA_records, poolA_arr)
		poolB_records = abind(poolB_records, poolB_arr)
	
		# Calculate convergence criteria if past burnin and at correct re-test interval
		if((step %% window_size == 0)&(step > burnin)){
			# Determine observations to use in calculation of R_hat
			window = floor(length((burnin+1):step)/2)
			use_obs = seq((step-window+1), step, thin)
			n_obs = length(use_obs)
		
			# Calculate community composition~env correlation for each chain
			corrstats = array(NA, dim=c(n_obs, nchains, 3, ncol(sites)), 
				dimnames=list(step=use_obs, chain=1:nchains, community=c('species','a','b'), env=1:ncol(sites)))
			for(n in use_obs){
			for(k in 1:nchains){
				corrstats[as.character(n), k, , ] = calc_envcorr(comm_records[,,k,n], topo_names, sites, c(), binary=F)
			}}

			# Calculate R_hat among chains
			Rhats = sapply(dimnames(corrstats)$community, function(y){
				sapply(1:ncol(sites), function(x){
					calc_Rhat(corrstats[,,y,x], alpha=.95)
				})
			})

			# Save records of Rhat
			if(exists('Rhat_records')){ 
				Rhat_records = abind(Rhat_records, Rhats)
			} else {
				Rhat_records = array(Rhats, dim=c(dim(Rhats),1), 
					dimnames=list(env=1:dim(Rhats)[1], community=dimnames(Rhats)[[2]]))
			}

			corr_flags = abs(Rhats-1) < Rhat_tol
		
			# Test whether all flags are TRUE
			converged = prod(corr_flags)>0
		}
		step = step + 1

	} # CLOSES while(converegence) loop
	} # CLOSES if(sim_mode=='converge')

	# Return simulation results in a list
	if(sim_mode=='converge') list(comm=comm_records, poolA=poolA_records, poolB=poolB_records, Tmat=Tmat_records, Rhat=Rhat_records)
	if(sim_mode=='fixed') list(comm=comm_records, poolA=poolA_records, poolB=poolB_records, Tmat=Tmat_records)
}


### Functions for checking status of simulation ###

# A function that calculates species abundances for all sites
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
#	type = 'species', 'a', 'b' : which partner should be tabulated	

calc_sa = function(comm, topo_names, type='species'){
	if(type=='species'){
		counts = apply(comm, 1, function(x) table(factor(x, 1:max(topo_names))))
	}
	
	if(type=='a'){
		counts = apply(comm, 1, function(x){
			table(factor(sapply(x, function(xi) ifelse(xi==0, 0, which(topo_names==xi, arr.ind=T)[1])), 1:nrow(topo_names)))
		})
	}

	if(type=='b'){
		counts = apply(comm, 1, function(x){
			table(factor(sapply(x, function(xi) ifelse(xi==0, 0, which(topo_names==xi, arr.ind=T)[2])), 1:ncol(topo_names)))
		})
	}

	t(counts)
}

# A function that calculates species richness for all sites
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
#	type = 'species', 'a', 'b' : which partner should be tabulated		
calc_rich = function(comm, topo_names, type){
	sa = calc_sa(comm, topo_names, type)
	apply(sa, 1, function(x) sum(x>0))
}

# A function that calculates species richness in the entire metacommunity 
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
#	type = 'species', 'a', 'b' : which partner should be tabulated		
calc_tot_rich = function(comm, topo_names, type){
	sa = calc_sa(comm, topo_names, type)
	sum(colSums(sa) > 0)
}

# A function that calculates the number of spaces occupied in each site
# 	comm = an N_C x N matrix of species present at each site
calc_occ = function(comm){
	apply(comm, 1, function(x) sum(x>0))
}

# A function that calculates the total number of spaces occupied in the entire metacommunity
# 	comm = an N_C x N matrix of species present at each site
calc_tot_occ = function(comm){
	sum(comm > 0)
}


# A function that calculates the compositional dissimilarity between all sites
# Standardization by total abundance is currently not implemented but could be easily added using vegan's decostand
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
#	type = 'species', 'a', 'b' : which partner should be tabulated
# 	metric = any dissimilarity method allowed by vegdist function in vegan 
#	binary = convert community matrix to presence/absence?	
calc_diss = function(comm, topo_names, type, metric, binary=F){
	comm_mat = calc_sa(comm, topo_names, type)
	if(binary) comm_mat = comm_mat > 0
	vegdist(comm_mat, method=metric)
}

# A function that calculates the correspondence between community and environmental variation
# Uses RDA on Hellinger-transformed species abundances
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
#	type = 'species', 'a', 'b' : which partner should be tabulated
#	env = matrix of environmental variables across sites (in same order as comm)
#	binary = convert community matrix to presence/absence?	
calc_rda = function(comm, topo_names, type, env, binary=F){
	# Calculate community matrix of species abudances
	comm_mat = calc_sa(comm, topo_names, type)

	# Convert to presence/absence if indicated
	if(binary) comm_mat = comm_mat > 0
	
	# Transform abundances
	comm_std = decostand(comm_mat, 'hellinger')

	# Conduct constrained ordination
	ord = rda(comm_std, env)

	# Return the R-square
	#RsquareAdj(ord)$adj.r.squared
	RsquareAdj(ord)$r.squared
}

# A function that calculates all community-environment correlation statistics
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
# 	metric = any dissimilarity method allowed by vegdist function in vegan 
#	binary = vector of logical indicated whether community matrix should be converted to presence/absence?
calc_envcorr = function(comm, topo_names, env, metric, binary){
	# Count the number of environmental axes
	Nenv = ifelse(is.null(dim(env)), 1, ncol(env))

	# Define array to hold statistics
	corr_mat = array(NA, dim=c(3, Nenv, length(metric)+3, length(binary)), 
		dimnames=list(community=c('species','a','b'), env=1:Nenv, measure=c('S','N','rda', metric), binary=binary))

	for(k in binary){
	for(type in c('species','a','b')){
	for(i in 1:Nenv){
		# Get environmental axis
		if(Nenv==1){ X = env } else { X = env[,i]}
		
		# Calculate correlation with species richness
		# NOTE: DOES NOT IMPLEMENT ABUNDANCE-WEIGHTED SPECIES DIVERSITY, SAME CALCULATION FOR BINARY=T/F
		corr_mat[type, i, 'S', k] = cor(X, calc_rich(comm, topo_names, type))

		# Calculate correlation with abundance
		# NOTE: THIS WILL BE THE SAME REGARDLESS OF TYPE OR BINARY
		corr_mat[type, i, 'N', k] = cor(X, calc_occ(comm))

		# Calculate RDA
		corr_mat[type, i, 'rda', as.character(k)] = calc_rda(comm, topo_names, type, X, k)
		
		# Calculate correlation for each dissimilarity metric
		for(j in metric){
			comm_diss = calc_diss(comm, topo_names, type, j, k)
			env_dist = dist(X)
			corr_mat[type, i, j, as.character(k)] = cor(comm_diss, env_dist, use='pairwise.complete.obs')
		}

	}}}
	
	corr_mat
}


# A function that calculates all community summary statistics
#	comm = an N_C x N matrix of species present at each site
#	topo_names = matrix describing mutualistic netiwork and numbering each partner
calc_commstats = function(comm, topo_names){

	# Calculate richness
	S_species = calc_rich(comm, topo_names, 'species')
	S_a = calc_rich(comm, topo_names, 'a')
	S_b = calc_rich(comm, topo_names, 'b')
	
	# Calculate total richness
	Stot_species = calc_tot_rich(comm, topo_names, 'species')
	Stot_a = calc_tot_rich(comm, topo_names, 'a')
	Stot_b = calc_tot_rich(comm, topo_names, 'b')

	# Calculate beta diversity
	Beta_species = Stot_species / mean(S_species, na.rm=T)
	Beta_a = Stot_a / mean(S_a, na.rm=T)
	Beta_b = Stot_b / mean(S_b, na.rm=T)

	# Calculate total abundance
	N = calc_occ(comm)
	Ntot = calc_tot_occ(comm)

	# Return dataframe
	# NOTE: tots will be the same across all communities
	data.frame(S_species, Stot_species, Beta_species, S_a, Stot_a, Beta_a, S_b, Stot_b, Beta_b, N, Ntot)
}


# A function that calculates interval-based scale reduction factor (from Brooks and Gelman 1998, DOI:10.1080/10618600.1998.10474787)
# THIS MAY NOT WORK WELL FOR DISCRETE QUANTITIES. MAY WANT TO DEFINE SOMETHING MORE EXACT BASED ON BINOMIAL DISTRIBUTION
#	x = matrix whose columns are independent chains and rows are realizations of the statistic to be compared
#	alpha = probability level for interval used to compare empirical distributions
calc_Rhat = function(x, alpha){
	within_chain = apply(x, 2, function(xi) diff(quantile(xi, c((1-alpha)/2, alpha + (1-alpha)/2))))
	all_chains = diff(quantile(x, c((1-alpha)/2, alpha + (1-alpha)/2)))
	R_hat = all_chains / mean(within_chain)

	R_hat
}





















