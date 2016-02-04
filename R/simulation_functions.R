## This script holds functions used by the community assembly of mutualists model (CAMM)
library(vegan) # vegdist

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
		add_a = sample(S_a, add_links)
		add_b = sample(S_b, add_links)
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
library(mvtnorm)
make_sites = function(N_C, rho_z){
	sd = c(1,1)
	mu = c(0,0)

	# Define covariance matrix
	S = matrix(c(sd[1],rho_z,rho_z,sd[2]), nrow=2, ncol=2)

	# Generate correlated gaussian rv using mvtnorm package
	sites = rmvnorm(mean=mu, sig=S, n=N_C)
	
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
#	S = number of species
#	nicheparms = a list of parameters controlling the shape of the distribution from which niches are sampled
#		mu = a vector of length 2 with the maximum niche optima
#		rho = the correlation between the niche optima
#		sigma = a vector of length 2 with the means of the gamma distributions
#		alpha = a vector of length 2 with the shape parameter of the gamma distributions
#		r = the correlation between niche widths

make_niches = function(N_S, nicheparms){
	
	# Generate niche optima as 2 correlated gaussian variables using mvtnorm packages
	S = matrix(c(1,nicheparms$rho,nicheparms$rho,1),2,2)
	mus_norm = rmvnorm(mean=c(0,0), sig=S, n=N_S)
	
	# Transform to uniform on (0,1)
	U = pnorm(mus_norm)

	# Transform to uniform on interval defined by nicheparms
	mu1 = qunif(U[,1], -nicheparms$mu[1], nicheparms$mu[1])
	mu2 = qunif(U[,2], -nicheparms$mu[2], nicheparms$mu[2])

	# Generate niche breadths as 2 correlated gaussian variables
	S = matrix(c(1, nicheparms$r, nicheparms$r, 1),2,2)
	sigmas_norm = rmvnorm(mean=c(0,0), sig=S, n=N_S)

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

# A function that generates a global speciesl abundance distribution.
# Species abundances may be correlated with niche optima or niche breadth or mutualist breadth.
# S = number of species
# distribution = list of parameters describing the functional form of the SAD
# condition = vector of quantities with which abundance should be correlated
# rho = strength of rank-order correlation between condition and abundance

make_gsad = function(S, distribution, condition=NA, rho=NA){
	require(poweRlaw)
	require(sads)
	require(permute)
	
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

	# THIS IS NOT VERY ROBUST AND SHOULD PROBABLY BE EDITED FOR BETTER ERROR HANDLING
	# MAY WANT TO RE-DRAW ABUNDANCES AFTER ALGORITHM FAILURE RATHER THAN IMMEDIATELY SENDING STOP ERROR MESSAGE
	if(!is.na(condition)){

		if(S != length(condition) ) stop('Number of species much match number of elements in condition')

		# The algorithm will swap the order of abundances until either
		#  the ranked correlation is within tol of the given rho
		# OR
		#  it has tried ntries swaps
		tol=0.01
		ntries = 10000

		c_rank = rank(condition)
		a_rank = rank(abuns)

		cor_func = function(x) cor(condition, abuns[x], method='spearman')
		use_order = 1:S
		this_rho = cor_func(use_order)
				
		n = 0
		while((abs(rho-this_rho) > tol) & (n < ntries)){
			n = n + 1
			swap = sample(S, 2, replace=F)
			new_order = use_order
			new_order[swap[1]] = use_order[swap[2]]
			new_order[swap[2]] = use_order[swap[1]]
			use_order = new_order
			this_rho = cor_func(use_order)
		}

		abuns = abuns[use_order]
		if(abs(rho-this_rho) > tol) stop(paste('Unable to correlate abundances and condition at rho=',rho,'\nConsider changing tolerance'))
	}
	
	abuns	
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

# A function that calculates the number of spaces occupied in each site
# 	comm = an N_C x N matrix of species present at each site
calc_occ = function(comm){
	apply(comm, 1, function(x) sum(x>0))
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
	corr_mat = array(NA, dim=c(3, Nenv, length(metric)+1, length(binary)), 
		dimnames=list(community=c('species','a','b'), env=1:Nenv, measure=c('rda', metric), binary=binary))

	for(k in binary){
	for(type in c('species','a','b')){
	for(i in 1:Nenv){
		# Get environmental axis
		if(Nenv==1){ X = env } else { X = env[,i]}

		# Calculate RDA
		corr_mat[type, i, 'rda', as.character(k)] = calc_rda(comm, topo_names, type, X, k)
		
		# Calculate correlation for each dissimilarity metric
		for(j in metric){
			comm_diss = calc_diss(comm, topo_names, type, j, k)
			env_dist = dist(X)
			corr_mat[type, i, j, as.character(k)] = cor(comm_diss, env_dist)
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

	# Calculate total abundance
	N = calc_occ(comm)

	# Return datafram
	data.frame(S_species, S_a, S_b, N)
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





















