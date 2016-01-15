## This script runs the community assembly of mutualists model (CAMM)

## Set options, load parameter values and simulation functions
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM'
setwd(working_dir)

options(stringsAsFactors=F)
source('parameter_file.R')
source('simulation_functions.R')

## Instantiate communities and mutualistic network
topo = make_topo(S_a, S_b, N_L, topology) # S_a x S_b matrix of association probabilities
topo_names = name_topo(topo) # S_a x S_b matrix with integers labeling specific associations
sites = make_sites(N_C, rho_z) # N_C x 2 matrix of env values
comm = make_comm(N_C, N, N_L) # initial N_C x N matrix of integers indicating which association is present
poolA = calc_pool(comm, topo_names, 1) # N_C x S_a matrix of presence of mutualist a in communities based on associations present
poolB = calc_pool(comm, topo_names, 2) # N_C x S_b matrix of presence of mutualist b in communities based on associations present

# Make random community for testing simulation
#comm = rand_comm(N_C, N, N_L)

## Generate random niches for each mutualist
niches_a = make_niches(S_a, nicheparms_a) # array of S_a x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
niches_b = make_niches(S_b, nicheparms_b)	# array of S_b x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)

# Plot niches for visualization
# plot_niches(niche_array, matrix of x-lims, site_env)
plot_niches(niches_a, matrix(c(-3,3,-3,3), 2, 2), add_env=sites)
plot_niches(niches_b, matrix(c(-3,3,-3,3), 2, 2), add_env=sites)


## A simulation of Markov process: dispersal implicit, rates independent and the same for both mutualists

# Generate empty array to hold changes in community and add initial community
comm_records = array(NA, c(N_C, N, reps+1)) 
comm_records[,,1] = comm
poolA_records = array(NA, c(N_C, S_a, reps+1)) 
poolA_records[,,1] = poolA
poolB_records = array(NA, c(N_C, S_b, reps+1))
poolB_records[,,1] = poolB
Tmat_records = array(NA, c(N_C, N_L+1, N_L+1, reps+1))

# Calculate N_C x N_L matrix of association probabilities at each site
assoc_probs = matrix(0, N_C, N_L)
for(i in 1:N_C){
for(j in 1:N_L){
	env = sites[i,]

	partners = which(topo_names==j, arr.ind=T)
	n_a = niches_a[partners[1],'mu',]
	n_b = niches_b[partners[2],'mu',]

	assoc_probs[i,j] = assoc_prob_func(env, n_a, n_b)
}}

x <- 1:n
PDF <- dls(x=x, N, alpha)
CDF <- pls(q=x, N=100, alpha=5)
par(mfrow=c(1,2))
plot(x,CDF, ylab="Cumulative Probability", type="b",
     main="Log-Series distribution, CDF")
plot(x,PDF, ylab="Probability", type="h",
     main="Log-Series distribution, PDF")
par(mfrow=c(1,1))

## Run simulation
for(step in 1:reps){
	# Mutualist mortality: mutualists present as associations cannot die
	poolA = die(comm, topo_names, poolA, mort_rate*mort_rate_a, 1)
	poolB = die(comm, topo_names, poolB, mort_rate*mort_rate_b, 2)

	# Mutualists disperse into communities independently and with probability based on niche-based filtering	
	poolA = disperse(sites, niches_a, poolA)
	poolB = disperse(sites, niches_b, poolB)

	# Calculate transition matrices for each site
	T_mat = calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate, assoc_probs) # array of square matrices of transition probabilities from one association to another
	Tmat_records[,,,step] = T_mat

	# Identify current state of each space and transition based on random numbers
	new_comm = matrix(NA, nrow=nrow(comm), ncol = ncol(comm))
	for(i in 1:nrow(comm)){
	for(j in 1:ncol(comm)){
		new_comm[i,j] = transition(T_mat[comm[i,j]+1,,i])
	}}
	comm = new_comm

	# Save communities
	comm_records[,,step+1] = comm
	poolA_records[,,step+1] = poolA
	poolB_records[,,step+1] = poolB

}

par(mfrow=c(5,1))
par(mar=c(4,4,1,1))
for(i in 1:5){
occupancy = colSums(comm_records[i,,]>0)
plot(occupancy, type='l', ylim=c(0,N))
abline(h=N, col=2)
}

par(mfrow=c(5,1))
par(mar=c(4,4,1,1))
for(i in 1:5){
richness = apply(comm_records[i,,], 2, function(x) length(unique(x)[unique(x)!=0]))
plot(richness, type='l', ylim=c(0,N_L))
abline(h=N_L, col=2)
}

# Presence of each partner
use_col = colorRampPalette(c('#90ABFC','#000000'))(S_a)
par(mfrow=c(5,1))
par(mar=c(4,4,1,1))
for(i in 1:5){
	plot(0,0, xlim=c(1,dim(comm_records)[3]), ylim=c(0,1), las=1, type='n', ylab=paste('Site',i), xlab='')
for(j in 1:S_a){
	lines(poolA_records[i,j,], col=use_col[j])
}
}

# Abundance spread of each species for a given site
use_col = colorRampPalette(c('#90ABFC','#000000'))(N_L)
par(mfrow=c(N_C,1))
for(i in 1:N_C){	
	plot(0,0, xlim=c(1,dim(comm_records)[3]), ylim=c(0,N), las=1, type='n', ylab=paste('Site',i), xlab='')
	counts = apply(comm_records[i,,], 2, function(x) table(factor(x, 1:10)))
	cum_counts = sapply(2:(N_L), function(i) colSums(counts[1:i,]))
	cum_counts = rbind(counts[1,], t(cum_counts))
	for(j in 1:S_a){
		lines(cum_counts[j,], col=use_col[j])
	}
}

## WORKING HERE
use_col = colorRampPalette(c('#90ABFC','#000000'))(S_b)
par(mfrow=c(N_C,1))
for(i in 1:N_C){	
	plot(0,0, xlim=c(1,dim(comm_records)[3]), ylim=c(0,N), las=1, type='n', ylab=paste('Site',i), xlab='')
	counts = apply(poolB_records[i,,], 2, function(x) table(factor(x, 1:10)))
	cum_counts = sapply(2:(S_b), function(i) colSums(counts[1:i,]))
	cum_counts = rbind(counts[1,], t(cum_counts))
	for(j in 1:S_b){
		lines(cum_counts[j,], col=use_col[j])
	}
}

use_col = colorRampPalette(c('#90ABFC','#000000'))(S_b)
par(mfrow=c(5,1))
par(mar=c(4,4,1,1))
for(i in 1:5){
	plot(0,0, xlim=c(1,dim(comm_records)[3]), ylim=c(0,1), las=1, type='n', ylab=paste('Site',i), xlab='')
for(j in 1:S_b){
	lines(poolB_records[i,j,], col=use_col[j])
}
}

cols = colorRampPalette(c('#90ABFC','#000000'))(5)
par(mfrow=c(1,2))
plot(0,0, type='n',, xlim=c(0,reps), ylim=c(0, S_a), ylab='S_a', xlab='Time') 
for(i in 1:5){
richness = colSums(poolA_records[i,,]>0)
points(richness, type='l', col=cols[i])
}
plot(0,0, type='n',, xlim=c(0,reps), ylim=c(0, S_b), ylab='S_b', xlab='Time') 
for(i in 1:5){
richness = colSums(poolB_records[i,,]>0)
points(richness, type='l', col=cols[i])
}




plot(comm_records[1,,]

