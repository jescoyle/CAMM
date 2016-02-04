## This script is used to analyze a simulation run

# Directory for saving figures
fig_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/Figures/'

# Load functions
code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'
source(paste(code_dir,'simulation_functions.R', sep=''))
source(paste(code_dir,'analysis_functions.R', sep=''))

# Load results
results_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/Runs/'
runID = 'testCONV2'

load(paste(results_dir, 'sim_results_',runID,'.RData', sep=''))
load(paste(results_dir, 'sim_object_',runID,'.RData', sep=''))
source(paste(results_dir, 'parameter_file_',runID,'.R', sep=''))

### Plot initial parameters

# Topology
svg(paste(fig_dir, runID, '_topology.svg', sep=''), height=4, width=4)
par(mar=c(0,0,0,0))
plot_topo(topo)
dev.off()

# Niches
svg(paste(fig_dir, runID, '_nichesA.svg', sep=''), height=4, width=8)
par(mar=c(4,4,0,0))
plot_niches(niches_a, grad=matrix(c(-3, 3, -3, 3), nrow=2), add_env=sites)
dev.off()

svg(paste(fig_dir, runID, '_nichesB.svg', sep=''), height=4, width=8)
par(mar=c(4,4,0,0))
plot_niches(niches_b, grad=matrix(c(-3, 3, -3, 3), nrow=2), add_env=sites)
dev.off()


### Code for simulations run for fixed number of time steps ###
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
par(mfrow=c(5,1))
for(i in 1:5){	
	plot(0,0, xlim=c(1,dim(comm_records)[3]), ylim=c(0,N), las=1, type='n', ylab=paste('Site',i), xlab='')
	counts = apply(comm_records[i,,], 2, function(x) table(factor(x, 1:N_L)))
	cum_counts = sapply(2:(N_L), function(j) colSums(counts[1:j,]))
	cum_counts = rbind(counts[1,], t(cum_counts))
	for(j in 1:S_a){
		lines(cum_counts[j,], col=use_col[j])
	}
}


use_col = colorRampPalette(c('#90ABFC','#000000'))(S_b)
par(mfrow=c(5,1))
for(i in 1:5){	
	plot(0,0, xlim=c(1,dim(comm_records)[3]), ylim=c(0,N), las=1, type='n', ylab=paste('Site',i), xlab='')
	
	counts = apply(comm_records[i,,], 2, function(x){
		table(factor(sapply(x, function(xi) ifelse(xi==0, 0, which(topo_names==xi, arr.ind=T)[2])), 1:S_b))
	})
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


### Code for simulations run until convergences (e.g. multiple chains) ###

## Compare simulation runs
thin=100

use_col = c('black','red','blue')
par(ask=T)


# Compare richness
use_obs = seq(1, dim(comm_records)[4], thin)

pdf(paste(fig_dir, runID, '_rich.pdf', sep=''), height=2*N_C, width=6)
layout(matrix(1:(2*N_C), ncol=2))
par(mar=c(2,4,1,1))
for(x in 1:2){
	site_order = order(sites[,x], decreasing=T)
for(y in c('a','b')){
for(i in site_order){
	plot(0,0, type='n', xlim=c(0, dim(comm_records)[4]), ylim=c(0,ifelse(y=='a',S_a,S_b)), xlab='', ylab=paste('Richness',y), las=1)
	richness = sapply(1:nchains, function(k) calc_rich(t(comm_records[i,,k,use_obs]), topo_names, y))
	for(k in 1:nchains){
		lines(use_obs, richness[,k], lwd=2, col=use_col[k])
	}
}
}}
dev.off()

# Compare richness ~ env correlations
# Note: we don't really expect correlations here.
thin=50
use_obs = seq(2,dim(comm_records)[4],thin)
use_col = c('black','red','blue')

pdf(paste(fig_dir, runID, '_rich_envcorr.pdf', sep=''), height=5, width=6)
par(mfrow=c(2,2))
par(mar=c(2,4,1,1))
for(x in 1:2){
for(y in c('a','b')){
	plot(0,0, type='n', xlim=c(0, dim(comm_records)[4]), ylim=c(-1,1), xlab='', ylab=paste('Richness',y,'~','Env',x), las=1)
	corr = sapply(1:nchains, function(k){
		sapply(use_obs, function(i){
			richness = calc_rich(comm_records[,,k,i], topo_names, y)
			cor(richness, sites[,x])
		})
	})
	for(k in 1:nchains) lines(use_obs, corr[,k], lwd=2, col=use_col[k])
}}
dev.off()

# Compare species composition ~ env correlations
thin=100
use_obs = seq(1,dim(comm_records)[4],thin)
use_col = c('black','red','blue')

pdf(paste(fig_dir, runID, '_envRDA.pdf', sep=''), height=5, width=6)
par(mfrow=c(2,2))
par(mar=c(2,4,1,1))
for(x in 1:2){
for(y in c('a','b')){
	plot(0,0, type='n', xlim=c(0, dim(comm_records)[4]), ylim=c(0,1), xlab='', ylab=paste('Mutualist',y,'~','Env',x), las=1)
	corr = sapply(1:nchains, function(k){
		sapply(use_obs, function(i){
			calc_rda(comm_records[,,k,i], topo_names, y, sites[,x], binary=F)
		})
	})
	for(k in 1:nchains) lines(use_obs, corr[,k], lwd=2, col=use_col[k])
}}
dev.off()

pdf(paste(fig_dir, runID, '_envRDA_presence.pdf', sep=''), height=5, width=6)
par(mfrow=c(2,2))
par(mar=c(2,4,1,1))
for(x in 1:2){
for(y in c('a','b')){
	plot(0,0, type='n', xlim=c(0, dim(comm_records)[4]), ylim=c(0,1), xlab='', ylab=paste('Mutualist',y,'~','Env',x), las=1)
	corr = sapply(1:nchains, function(k){
		sapply(use_obs, function(i){
			calc_rda(comm_records[,,k,i], topo_names, y, sites[,x], binary=T)
		})
	})
	for(k in 1:nchains) lines(use_obs, corr[,k], lwd=2, col=use_col[k])
}}
dev.off()


# Plot Rhats
rhat_max = max(Rhat_records)
nrecs = dim(comm_records)[4]

pdf(paste(fig_dir, runID, '_Rhat_envRDA.pdf', sep=''), height=4, width=5)
par(mfrow=c(2,2))
par(mar=c(2,4,1,1))
for(x in 1:2){
for(y in c('a','b')){
	plot(0,0, type='n', las=1, xlim=c(0,nrecs), ylim=c(0, rhat_max), xlab='', ylab=paste('Mutualist',y,'~ Env', x))
	yvals = Rhat_records[x,y,]
	xvals = 150 + (0:(length(yvals)-1))*50
	lines(xvals, yvals, lwd=2)
	abline(h=1)
}}
dev.off()

