## This script contains functions for visualizing, summarizing and analyzing simulation results

## CURRENTLY CODE IS NOT WRITTEN AS FUNCTIONS
 
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


use_col = colorRampPalette(c('#90ABFC','#000000'))(S_b)
par(mfrow=c(N_C,1))
for(i in 1:N_C){	
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

