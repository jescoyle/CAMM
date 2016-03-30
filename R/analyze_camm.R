## This script is used to analyze a simulation runs

options(stringsAsFactors=F)
library(reshape2)
library(lattice) # Plots

# Main directory
working_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/'
setwd(working_dir)

# Directory for saving figures
fig_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/Figures/'

# Location of results
results_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/Runs/Summaries_5/'

# Load functions
code_dir = 'C:/Users/jrcoyle/Documents/UNC/Projects/CAMM/GitHub/R/'
#source(paste(code_dir,'simulation_functions.R', sep=''))
source(paste(code_dir,'analysis_functions.R', sep=''))

## Check which results are done
#done = read.table('obs_done.txt', header=F)
#
#parms_done = unique(gsub('_[0-9]+_results.RData', '', done$V1))
#parms_done = unique(gsub('_[0-9]+_simobject.RData', '', done$V1))
#
#parms_done_df = data.frame(t(sapply(parms_done, get_parms)))
#parms_done_df$runs_done = sapply(parms_done, function(x) length(grep(paste0(x, '_[0-9]+_results.RData'), done$V1)))

### Analyze multiple runs across a set of parameters ###

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

	# Bind values to data
	this_data = cbind(parm_vals, this_data)

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

	# Bind values to data
	this_data = cbind(parm_vals, this_data)

	# Add to growing data frame
	cor_summary = rbind(cor_summary, this_data)
}


## Run 1: incrementing over stregth of mutualism (omega = o) and relative mortality of unassociated mutualists (mort_rate_a = mra, mort_rate_b = mrb)

cor_summary$o = as.numeric(cor_summary$o)
cor_summary$mra = as.numeric(cor_summary$mra)
cor_summary$mrb = as.numeric(cor_summary$mrb)

## Plot correlations between mutualist communities and environment
## error bars show 95th percentile from 100 different starts

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
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()

# Results from obligate mutualism only (omega=1)
plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & o==1)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
jit_fact = .1
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'obligate mutualism_mort_rates_RDAmean.pdf'), height=3, width=9)
xyplot(cor_a ~ factor(mrb) | mra, groups = env, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab='Relative mortality of unassociated mutualist B', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		y = as.numeric(y)
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a]), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()

# RDA within chain variance (within 10 chains)
plot_data = subset(cor_summary, measure=='rda' & summary=='var')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
jit_fact = 0.008
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'mutualism_strength_vs_mort_rates_RDAvar.pdf'), height=7, width=9)
xyplot(cor_a ~ o | mrb + mra, groups = env, data=means, ylim=c(-0.001, .02),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
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
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
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


## IT MIGHT ALSO BE USEFUL TO EXAMINE THE SYMBIONT POOLS SEPARATELY
## This section run in interactive mode on killdevil cluster to avoid transfering results files
## SECTION NOT DONE OR COMPLETE

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



## Run 2: incrementing over strength of mutualism (omega = o), topology (topo) and type of environmental filtering (envfilt)


## Plot correlations between mutualist communities and environment
## error bars show 95th percentile from 100 different starts

cor_summary$o = as.numeric(cor_summary$o)
cor_summary$topo = factor(cor_summary$topo, levels=c('one2one','one2many','many2many'))
cor_summary$envfilt = factor(cor_summary$envfilt, levels=c('none','same','opposite','all'))

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

pdf(paste0(fig_dir, 'mutualism_strength_vs_topo&envfilt_RDAmean.pdf'), height=7, width=9)
xyplot(cor_a ~ o | topo + envfilt, groups = env, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=c('topology','env filters'), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()

# Plot for only obligate mutualism
plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & o==1)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
jit_fact = 0.1
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'topo_vs_envfilt_RDAmean_o=1.pdf'), height=5, width=7.5)
xyplot(cor_a ~ topo | envfilt, groups = env, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab='Network Topology', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(as.numeric(x)+jit_a[groups[subscripts]], low95s$cor_a[subscripts], as.numeric(x)+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(as.numeric(x)+jit_b[groups[subscripts]], low95s$cor_b[subscripts], as.numeric(x)+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(as.numeric(x)+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(as.numeric(x)+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=c('env filters'), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
dev.off()


### Run 3: Incrementng over topology, direction and strength of env filtering and 

# sigA : std dev of niche breadth for A
# sigB : std dev of niche breadth for B
# topo : type of topology
# envfilt : none, opposite (A ~ 2, B ~ 1), same (A ~ 2, B ~ 2), all

cor_summary$sigA = as.numeric(cor_summary$sigA)
cor_summary$sigB = as.numeric(cor_summary$sigB)
cor_summary$topo = factor(cor_summary$topo, levels = c('one2one','one2many','many2many'))
cor_summary$envfilt = factor(cor_summary$envfilt, levels = c('none','same','opposite','all'))
comm_summary$sigA = as.numeric(comm_summary$sigA)
comm_summary$sigB = as.numeric(comm_summary$sigB)
comm_summary$topo = factor(comm_summary$topo, levels = c('one2one','one2many','many2many'))
comm_summary$envfilt = factor(comm_summary$envfilt, levels = c('none','same','opposite','all'))


a_pch = c(16, 1)
b_pch = c(15, 0)
jit_fact = 0.15
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

# RDA within chain mean

pdf(paste0(fig_dir, 'envfilt_strength_by_type_RDAmean.pdf'), height=8.5, width=11)

for(t in levels(cor_summary$topo)){

plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & topo==t)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

print(
xyplot(cor_a ~ factor(sigA) | sigB + envfilt, groups = env, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(RDA~~R^2),main=paste('Topology =',t),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
)
}
dev.off()

plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & env==2)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
use_col = c('black','red','cornflowerblue')
pch_1 = c(0,1,2,5)

pdf(paste0(fig_dir, 'cora_vs_corb_mean_across_envfilt2_type_strength.pdf'), height=10, width=10)
xyplot(cor_a  ~ cor_b | sigA + sigB, groups=envfilt, data=means,
	ylim=c(0,.65), xlim=c(0,.65), ylab=expression(S[a]), xlab=expression(S[b]),
	panel=function(x, y, subscripts, groups){
		panel.segments(x, low95s$cor_a[subscripts], x, up95s$cor_a[subscripts], col='grey')
		panel.segments(low95s$cor_b[subscripts], y, up95s$cor_b[subscripts], y, col='grey')
		colvec = use_col[as.numeric(means$topo[subscripts])]
		panel.xyplot(x, y, pch=pch_1[groups[subscripts]], col=colvec)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[a], sigma[b]), fg='transparent'),
	key=list(space='right', points=list(pch=pch_1), text=list(levels(means$envfilt)))
)
dev.off()


pdf(paste0(fig_dir, 'envfilt_strength_by_type_Smean_cor.pdf'), height=8.5, width=11)

for(t in levels(cor_summary$topo)){

plot_data = subset(cor_summary, measure=='S' & summary=='mean' & topo==t)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')


print(
xyplot(cor_a ~ factor(sigA) | sigB + envfilt, groups = env, data=means, ylim = c(-1,1),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(RDA~~R^2),main=paste('Topology =',t),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
)
}
dev.off()

# Diversity and abundance changes?

plot_data = subset(comm_summary, summary=='mean')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
jit_fact = 0.18
jit = c(-1, 0, 1)*jit_fact

use_pch = c(0,1,2)

# Might want to plot these as percentages of total possible S
pdf(paste0(fig_dir, 'envfilt_strength_by_type_Sa_mean.pdf'), height=8.5, width=11)
xyplot(S_a ~ factor(sigA) | sigB + envfilt, groups=topo, data=means, ylim=c(0,32),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(Mean~~S[A]),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit, low95s$S_a[subscripts], x+jit, up95s$S_a[subscripts])
		panel.xyplot(x+jit, y, pch=use_pch, col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()
pdf(paste0(fig_dir, 'envfilt_strength_by_type_Sb_mean.pdf'), height=8.5, width=11)
xyplot(S_b ~ factor(sigA) | sigB + envfilt, groups=topo, data=means, ylim=c(0,32),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(Mean~~S[B]),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit, low95s$S_b[subscripts], x+jit, up95s$S_b[subscripts])
		panel.xyplot(x+jit, y, pch=use_pch, col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()
pdf(paste0(fig_dir, 'envfilt_strength_by_type_N_mean.pdf'), height=8.5, width=11)
xyplot(N ~ factor(sigA) | sigB + envfilt, groups=topo, data=means, ylim=c(0,110),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab='Mean N',
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit, low95s$N[subscripts], x+jit, up95s$N[subscripts])
		panel.xyplot(x+jit, y, pch=use_pch, col=1)
		panel.abline(h=100, lty=3)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()


## Sa vs Sb vs RDA R2
plot_data = subset(comm_summary, summary=='mean'&topo!='one2one')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
pch_1 = c(0,1,2,5)
use_col = c(1:2)

pdf(paste0(fig_dir, 'Sa_vs_Sb_mean_across_envfilt_type_strength.pdf'), height=10, width=10)
xyplot(S_a  ~ S_b | sigA + sigB, groups=envfilt, data=means,
	ylim=c(0,32), xlim=c(0,12), ylab=expression(S[a]), xlab=expression(S[b]),
	panel=function(x, y, subscripts, groups){
		panel.segments(x, low95s$S_a[subscripts], x, up95s$S_a[subscripts], col='grey')
		panel.segments(low95s$S_b[subscripts], y, up95s$S_b[subscripts], y, col='grey')
		panel.xyplot(x, y, pch=pch_1[groups[subscripts]], col=1+as.numeric(means$topo[subscripts]=='many2many'))
		panel.abline(h=30, lty=3)
		panel.abline(v=10, lty=3)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[a], sigma[b]), fg='transparent'),
	key=list(space='right', points=list(pch=pch_1), text=list(levels(means$envfilt)))

)
dev.off()

plot_data = subset(cor_summary, measure=='rda')
plot_data = merge(plot_data, comm_summary)
plot_data = subset(plot_data, summary=='mean'&env==1)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

plot(cor_b ~ S_b, data=means)




### RUN 4: Run out to 8000 time steps

cor_summary$T = as.numeric(substring(cor_summary$time, 2))
comm_summary$T = as.numeric(substring(comm_summary$time, 2))
cor_summary$topo = factor(cor_summary$topo, levels = c('one2one','one2many','many2many'))
cor_summary$envfilt = factor(cor_summary$envfilt, levels = c('none','same','opposite','all'))
comm_summary$topo = factor(comm_summary$topo, levels = c('one2one','one2many','many2many'))
comm_summary$envfilt = factor(comm_summary$envfilt, levels = c('none','same','opposite','all'))


## Plot changes through time for each topology and envfilt
plot_data = subset(cor_summary, summary=='mean'&measure=='rda')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

a_pch = c(16, 1)
b_pch = c(15, 0)
jit_fact = 100
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'envfilt_by_topo_through_time_RDAmean.pdf'), height=8.5, width=11)
print(
xyplot(cor_a ~ T | topo + envfilt, groups = env, data=means, ylim = c(0,0.6),
	scales=list(alternating=1), xlab='Time (# steps)', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x[groups==1]+jit_a[1], y[groups==1], type='l', col='grey')
		panel.xyplot(x[groups==2]+jit_a[2], y[groups==2], type='l', col='grey')
		panel.xyplot(x[groups==1]+jit_b[1], means$cor_b[subscripts][groups==1], col='grey', type='l')
		panel.xyplot(x[groups==2]+jit_b[2], means$cor_b[subscripts][groups==2], col='grey', type='l')
		
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)

	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=c('Topology','Filtering'), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
)
dev.off()

plot_data = subset(cor_summary, summary=='var'&measure=='rda')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'envfilt_by_topo_through_time_RDAvar.pdf'), height=8.5, width=11)
print(
xyplot(cor_a ~ T | topo + envfilt, groups = env, data=means, ylim=c(0, max(up95s[,c('cor_a','cor_b')])),
	scales='free', xlab='Time (# steps)', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x[groups==1]+jit_a[1], y[groups==1], type='l', col='grey')
		panel.xyplot(x[groups==2]+jit_a[2], y[groups==2], type='l', col='grey')
		panel.xyplot(x[groups==1]+jit_b[1], means$cor_b[subscripts][groups==1], col='grey', type='l')
		panel.xyplot(x[groups==2]+jit_b[2], means$cor_b[subscripts][groups==2], col='grey', type='l')
	
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)	
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=c('Topology','Filtering'), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
)
dev.off()

plot_data = subset(comm_summary, summary=='mean')

# Scale richness of symbiont by total possible species
plot_data[plot_data$topo=='one2one',c('S_b', 'Stot_b')] = plot_data[plot_data$topo=='one2one',c('S_b', 'Stot_b')] / 30
plot_data[plot_data$topo!='one2one',c('S_b', 'Stot_b')] = plot_data[plot_data$topo!='one2one',c('S_b', 'Stot_b')] / 10

means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
jit_fact = 100
jit = c(-1, 0, 1)*jit_fact

use_pch = c(0,1,2)
names(use_pch) = levels(plot_data$topo)


pdf(paste0(fig_dir, 'envfilt_by_topo_through_time_Sa_mean.pdf'), height=5, width=7)
xyplot(S_a ~ T | envfilt, groups=topo, data=means, ylim=c(0,32), xlim=c(0,9000),
	scales=list(alternating=1), xlab='Time (#steps)', ylab=expression(Mean~~S[A]),
	panel=function(x, y, subscripts, groups){
		for(i in 1:3){
			i_name = levels(plot_data$topo)[i]
			panel.xyplot(x[groups[subscripts]==i_name]+jit[i], y[groups[subscripts]==i_name], type='l', col='grey')
		}
		panel.segments(x+jit[groups[subscripts]], low95s$S_a[subscripts], x+jit[groups[subscripts]], up95s$S_a[subscripts])
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name='Filtering', fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()

pdf(paste0(fig_dir, 'envfilt_by_topo_through_time_Sb_mean.pdf'), height=5, width=7)
xyplot(S_b ~ T | envfilt, groups=topo, data=means, ylim=c(0,1.1), xlim=c(0,9000),
	scales=list(alternating=1), xlab='Time (#steps)', ylab=expression(Mean~~S[B]~~Prop.),
	panel=function(x, y, subscripts, groups){
		for(i in 1:3){
			i_name = levels(plot_data$topo)[i]
			panel.xyplot(x[groups[subscripts]==i_name]+jit[i], y[groups[subscripts]==i_name], type='l', col='grey')
		}
		panel.segments(x+jit[groups[subscripts]], low95s$S_b[subscripts], x+jit[groups[subscripts]], up95s$S_b[subscripts])
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name='Filtering', fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()

pdf(paste0(fig_dir, 'envfilt_by_topo_through_time_N_mean.pdf'), height=5, width=7)
xyplot(N ~ T | envfilt, groups=topo, data=means, ylim=c(0,101), xlim=c(0,9000),
	scales=list(alternating=1), xlab='Time (#steps)', ylab='Mean N',
	panel=function(x, y, subscripts, groups){
		for(i in 1:3){
			i_name = levels(plot_data$topo)[i]
			panel.xyplot(x[groups[subscripts]==i_name]+jit[i], y[groups[subscripts]==i_name], type='l', col='grey')
		}
		panel.segments(x+jit[groups[subscripts]], low95s$N[subscripts], x+jit[groups[subscripts]], up95s$N[subscripts])
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name='Filtering', fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()


### RUN 5: Effect of changing number of symbiont and host species
cor_summary$topo = factor(cor_summary$topo, levels = c('one2many','many2many'))
cor_summary$envfilt = factor(cor_summary$envfilt, levels = c('same','opposite'))
comm_summary$topo = factor(comm_summary$topo, levels = c('one2many','many2many'))
comm_summary$envfilt = factor(comm_summary$envfilt, levels = c('same','opposite'))
cor_summary$a2b = cor_summary$Sa / cor_summary$Sb
cor_summary$NL2a = sapply(cor_summary$NL / cor_summary$Sa, match_NLratio)
comm_summary$a2b = comm_summary$Sa / comm_summary$Sb
comm_summary$NL2a = sapply(comm_summary$NL / comm_summary$Sa, match_NLratio)

for(i in c('Sa','Sb','NL')){
	cor_summary[,i] = as.numeric(cor_summary[,i])
	comm_summary[,i] = as.numeric(comm_summary[,i])
}

NLratios = c(1.25,1.5,2,3)
match_NLratio = function(x){
	ds = abs(x - NLratios)
	NLratios[ds==min(ds)]
}

plot_data = subset(cor_summary, topo=='one2many'&summary=='mean'&measure=='rda')
plot_data$a2b = plot_data$Sa / plot_data$Sb
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

a_pch = c(16, 1)
b_pch = c(15, 0)
jit_fact = 2
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)


pdf(paste0(fig_dir, 'Sratios_topo-one2many_RDAmean.pdf'), height=5, width=10)
print(
xyplot(cor_a ~ Sa | a2b + envfilt, groups = env, data=means, ylim = c(0,1),
	scales=list(alternating=1), xlab=expression(S[A]), ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)

	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(S[A]:S[B],Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
		text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
)
)
dev.off()


jit = c(-.6,.6)

pdf(paste0(fig_dir, 'Sratios_topo-many2many_RDAmean.pdf'), height=7, width=12)

for(v in unique(cor_summary$envfilt)){

plot_data = subset(cor_summary, topo=='many2many'&summary=='mean'&measure=='rda'&envfilt==v)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

for(yvar in c('cor_a','cor_b')){
if(yvar=='cor_a'){
	use_pch = a_pch
	partner = 'Host'
}
if(yvar=='cor_b'){
	use_pch = b_pch
	partner = 'Symbiont'
}

print(
xyplot(means[,yvar] ~ Sb | a2b + NL2a, groups = env, data=means, ylim = c(0,1), main=paste('Filtering =',v,' : ',partner),
	scales=list(alternating=1), xlab=expression(S[B]), ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit[groups[subscripts]], low95s[subscripts,yvar], x+jit[groups[subscripts]], up95s[subscripts,yvar])
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(S[A]:S[B],N[L]:S[A]), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(expression(" ~ "*E[1]," ~ "*E[2])))
)
)

}}
dev.off()

# Hold Sb constant
pdf(paste0(fig_dir, 'Sratios_topo-many2many_RDAmean_compareSa.pdf'), height=9, width=6)

for(b in unique(cor_summary$Sb)){

plot_data = subset(cor_summary, topo=='many2many'&summary=='mean'&measure=='rda'&Sb==b)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

for(yvar in c('cor_a','cor_b')){
if(yvar=='cor_a'){
	use_pch = a_pch
	partner = 'Host'
}
if(yvar=='cor_b'){
	use_pch = b_pch
	partner = 'Symbiont'
}
print(
xyplot(means[,yvar] ~ Sa | envfilt + NL2a, groups = env, data=means, ylim = c(0,1), main=bquote(.(partner)~~(S[B]==.(b))),
	scales=list(alternating=1), xlab=expression(S[A]), ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit[groups[subscripts]], low95s[subscripts,yvar], x+jit[groups[subscripts]], up95s[subscripts,yvar])
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(Filtering, N[L]:S[A]), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(expression(" ~ "*E[1]," ~ "*E[2])))
)
)

}}
dev.off()

jit = c(-.01,.01)
pdf(paste0(fig_dir, 'Sratios_topo-many2many_RDAmean_compareNL.pdf'), height=9, width=6)

for(b in unique(cor_summary$Sb)){

plot_data = subset(cor_summary, topo=='many2many'&summary=='mean'&measure=='rda'&Sb==b)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

for(yvar in c('cor_a','cor_b')){
if(yvar=='cor_a'){
	use_pch = a_pch
	partner = 'Host'
}
if(yvar=='cor_b'){
	use_pch = b_pch
	partner = 'Symbiont'
}
print(
xyplot(means[,yvar] ~ NL2a | envfilt + Sa, groups = env, data=means, ylim = c(0,1), main=bquote(.(partner)~~(S[B]==.(b))),
	scales=list(alternating=1), xlab=expression(N[L]:S[A]), ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		panel.segments(x+jit[groups[subscripts]], low95s[subscripts,yvar], x+jit[groups[subscripts]], up95s[subscripts,yvar])
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(Filtering, S[A]), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(expression(" ~ "*E[1]," ~ "*E[2])))
)
)

}}
dev.off()




###########################################################################
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

