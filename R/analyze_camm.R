## This script is used to analyze a simulation runs

options(stringsAsFactors=F)
library(reshape2)

# Directory for saving figures
fig_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/Figures/'

# Location of results
results_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/Runs/Summaries_1/'

# Load functions
code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'
#source(paste(code_dir,'simulation_functions.R', sep=''))
source(paste(code_dir,'analysis_functions.R', sep=''))

## Check which results are done
done = read.table('done.txt', header=F)

parms_done = unique(gsub('_[0-9]+_results.RData', '', done$V1))
parms_done_df = data.frame(t(sapply(parms_done, get_parms)))
parms_done_df$runs_done = sapply(parms_done, function(x) length(grep(paste0(x, '_[0-9]+_results.RData'), done$V1)))

### Analyze multiple runs across a set of parameters ###

## Run incrementing over stregth of mutualism (omega = o) and relative mortality of unassociated mutualists (mort_rate_a = mra, mort_rate_b = mrb)

# Get list of community richness and abundance summaries
comm_filelist = list.files(results_dir, 'comm_summary.csv')

# Go through files and append together into a data frame
comm_summary = data.frame()
for(f in comm_filelist){
	this_data = read.csv(paste0(results_dir, f))
	
	# Extract run name
	runID = gsub('_comm_summary.csv', '', f)

	# Extract parameter values
	parm_vals = get_parms(runID)

	# Load results


	# Bind values to data
	this_data = cbind(t(parm_vals), this_data)

	# Add to growing data frame
	comm_summary = rbind(comm_summary, this_data)
}

# Get list of community richness and abundance summaries
cor_filelist = list.files(results_dir, 'cor_summary.csv')

# Go through files and append together into a data frame
cor_summary = data.frame()
for(f in cor_filelist){
	this_data = read.csv(paste0(results_dir, f))
	
	# Extract run name
	runID = gsub('_cor_summary.csv', '', f)

	# Extract parameter values
	parm_vals = get_parms(runID)

	# Get environmental values for each site


	# Bind values to data
	this_data = cbind(t(parm_vals), this_data)

	# Add to growing data frame
	cor_summary = rbind(cor_summary, this_data)
}

## Plot correlations between mutualist communities and environment
## error bars show 95th percentile from 100 different starts
library(lattice)

a_pch = c(16, 1)
b_pch = c(15, 0)
jit_fact = 0.008
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

# RDA within chain mean
plot_data = subset(cor_summary, measure=='rda' & summary=='mean')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'mutualism_strength_vs_mort_rates_RDAmean.pdf'), height=7, width=9)
xyplot(cor_a ~ o | mrb + mra, groups = env, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups], low95s$cor_a[subscripts], x+jit_a[groups], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups], low95s$cor_b[subscripts], x+jit_b[groups], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups], y, pch=a_pch[groups], col=1)
		panel.xyplot(x+jit_b[groups], means$cor_b[subscripts], pch=b_pch[groups], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()

# RDA within chain variance (within 10 chains)
plot_data = subset(cor_summary, measure=='rda' & summary=='var')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'mutualism_strength_vs_mort_rates_RDAvar.pdf'), height=7, width=9)
xyplot(cor_a ~ o | mrb + mra, groups = env, data=means, ylim=c(-0.001, .02),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups], low95s$cor_a[subscripts], x+jit_a[groups], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups], low95s$cor_b[subscripts], x+jit_b[groups], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups], y, pch=a_pch[groups], col=1)
		panel.xyplot(x+jit_b[groups], means$cor_b[subscripts], pch=b_pch[groups], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()

# Jaccard within chain mean
plot_data = subset(cor_summary, measure=='jaccard' & summary=='mean')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'mutualism_strength_vs_mort_rates_Jaccardmean.pdf'), height=7, width=9)
xyplot(cor_a ~ o | mrb + mra, groups = env, data=means, ylim=c(-1,1),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab='Jaccard Similarity ~ Env',
	panel=function(x, y, subscripts, groups){
		panel.abline(h=0, col='grey')
		panel.segments(x+jit_a[groups], low95s$cor_a[subscripts], x+jit_a[groups], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups], low95s$cor_b[subscripts], x+jit_b[groups], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups], y, pch=a_pch[groups], col=1)
		panel.xyplot(x+jit_b[groups], means$cor_b[subscripts], pch=b_pch[groups], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()

# Species richness and abundance
# ACTUALLY THIS DOESN'T MAKE SENSE B/C COMMUNITIES ARE NOT EXPECTED TO HAVE THE SAME ENV ACROSS SIMULATIONS
plot_data = subset(comm_summary, summary=='mean')
plot_data = plot_data[order(plot_data$o),]
means = subset(plot_data, stat=='mean')
use_col = rainbow(20)

xyplot(S_a ~ o | mrb + mra, groups = comm, data=means, type='l',ylim=c(0,30),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(S[A]),
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent')
)

xyplot(S_b ~ o | mrb + mra, groups = comm, data=means, type='l',ylim=c(0,10),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(S[A]),
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent')
)

xyplot(N ~ o | mrb + mra, groups = comm, data=means, type='l', ylim=c(0,100),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(S[A]),
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent')
)


## IT WOULD ALSO BE USEFUL TO EXAMINE THE SYMBIONT POOLS SEPARATELY
## This section run in interactive mode on killdevil cluster to avoid transfering results files

results_dir = './Results/'
source('analysis_functions.R')


# Get list of community richness and abundance summaries
results_filelist = list.files(results_dir, 'results.RData')
runs = unique(gsub('_[0-9]+_results.RData', '', results_filelist))


# Go through files and append together into a data frame
for(this_run in runs){
	these_files = results_filelist[grep(paste0(this_run,'_[0-9]+'), results_filelist)]

	poolA_summary = data.frame()
	poolB_summary = data.frame()
	
	for(f in these_files){
		load(paste0(results_dir, f))

		poolA = end_metacomms['poolA',]
		poolA = acast(melt(poolA), Var1 ~ Var2 ~ L1)

		poolB = end_metacomms['poolB',]
		poolB = acast(melt(poolB), Var1 ~ Var2 ~ L1)

		A_probs = apply(poolA, c(1,2), function(x) sum(x)/nchains)
		B_probs = apply(poolB, c(1,2), function(x) sum(x)/nchains)

		A_probs = melt(A_probs)
		B_probs = melt(B_probs)

		# Extract run name
		runID = gsub('_results.RData', '', f)

		# Extract parameter values
		parm_vals = get_parms(runID)

		# Bind values to data
		A_data = cbind(t(parm_vals[1:3]), names(parm_vals)[4], A_probs)
		names(A_data)[4:7] = c('run','comm','species','prob')

		B_data = cbind(t(parm_vals[1:3]), names(parm_vals)[4], B_probs)
		names(B_data)[4:7] = c('run','comm','species','prob')

		# Add to growing data frame
		poolA_summary = rbind(poolA_summary, A_data)
		poolB_summary = rbind(poolB_summary, B_data)

	}

	write.csv(poolA_summary, paste0(results_dir, this_run, '_poolA.csv'), row.names=F)
	write.csv(poolB_summary, paste0(results_dir, this_run, '_poolB.csv'), row.names=F)

	poolA_arr = acast(poolA_summary, comm~species~run, value.var='prob') 
	poolB_arr = acast(poolB_summary, comm~species~run, value.var='prob')

	poolA_stats = apply(poolA_arr, c(1,2), function(x) c(mean=mean(x, na.mrm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975))))
	poolB_stats = apply(poolB_arr, c(1,2), function(x) c(mean=mean(x, na.mrm=T), var=var(x, na.rm=T), quantile(x, c(0.025, 0.5, 0.975))))

	# Calculate correlations and richness



### Analyze a single run ###

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

