## This script is used to analyze a simulation runs

options(stringsAsFactors=F)
library(reshape2)
library(lattice) # Plots
library(latticeExtra)

# Main directory
working_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/'
setwd(working_dir)

# Directory for saving figures
fig_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/Figures/'

# Location of results
results_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/Runs/Summaries_1-new/'

# Load functions
code_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/GitHub/CAMM/R/'
#source(paste(code_dir,'simulation_functions.R', sep=''))
source(paste(code_dir,'analysis_functions.R', sep=''))
source(paste(code_dir,'simulation_functions.R', sep=''))

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


## In each run, evaluate effects of parameters on:
## 1) Mean richness of hosts and symbionts
## 2) Turnover of hosts and symbionts (total richness / mean richness)
## 3) Correlation between host and symbiont richness
## 4) Mean abundance
## 5) Correlation between environment and host & symbiont community structure
## 6) Correlation between environment and host & symbiont richness
## 7) Correlation between environment and abundance 

## Define symbols used for plotting throughout
a_pch = c(16, 1)
b_pch = c(15, 0)
n_pch = c(17, 2)


## Run 1: incrementing over stregth of mutualism (omega = o) and relative mortality of unassociated mutualists (mort_rate_a = mra, mort_rate_b = mrb)

cor_summary$o = as.numeric(cor_summary$o)
cor_summary$mra = as.numeric(cor_summary$mra)
cor_summary$mrb = as.numeric(cor_summary$mrb)
comm_summary$o = as.numeric(comm_summary$o)
comm_summary$mra = as.numeric(comm_summary$mra)
comm_summary$mrb = as.numeric(comm_summary$mrb)


## 1) Mean richness of hosts and symbionts

plot_data = subset(comm_summary, summary=='mean')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

jit = 0.01*c(-1,1)

pdf(paste0(fig_dir, 'RUN 1/', 'mutualism_strength_vs_mort_rates_mean_S.pdf'), height=7, width=9)
lp1 = xyplot(S_a ~ o | mrb + mra, data=means, ylim = c(-.1,30),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(Mean~~S[A]),
	panel=function(x, y, subscripts){
		panel.segments(x+jit[1], low95s$S_a[subscripts], x+jit[1], up95s$S_a[subscripts], col=1)
		panel.xyplot(x+jit[1], y, pch=a_pch[1], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent')
)

lp2 = xyplot(S_b ~ o | mrb + mra, data=means, ylim = c(-.1,10),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab=expression(Mean~~S[B]),
	panel=function(x, y, subscripts){
		panel.segments(x, low95s$S_b[subscripts], x, up95s$S_b[subscripts], col=1)
		panel.xyplot(x, y, pch=b_pch[1], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent')
)
doubleYScale(lp1, lp2, add.axis=T, add.ylab2=T, text=expression(S[A], S[B]))
dev.off()

## 2) Turnover of hosts and symbionts (total richness / mean richness)

## 3) Correlation between host and symbiont richness

## 4) Mean abundance
plot_data = subset(comm_summary, summary=='mean')
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'RUN 1/', 'mutualism_strength_vs_mort_rates_mean_N.pdf'), height=7, width=9)
xyplot(N ~ o | mrb + mra, data=means, ylim=c(80,100),
	scales=list(alternating=1), xlab='Strength of mutalism (omega)', ylab='Mean Total Abundance',
	panel=function(x, y, subscripts){
		panel.segments(x, low95s$N[subscripts], x, up95s$N[subscripts], col=1)
		panel.xyplot(x, y, pch=n_pch[1], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a],m[b]), fg='transparent')
)
dev.off()

# From obligate mutualism onlye
plot_data = subset(comm_summary, summary=='mean' & o==1)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'RUN 1/', 'mutualism_strength_vs_mort_rates_mean_N_o=1.pdf'), height=3, width=4)
xyplot(N ~ mrb | mra, data=means, ylim=c(80,100), layout=c(3,1),
	scales=list(alternating=1), xlab=expression(m[b]), ylab='Mean Total Abundance',
	panel=function(x, y, subscripts){
		panel.segments(x, low95s$N[subscripts], x, up95s$N[subscripts], col=1)
		panel.xyplot(x, y, pch=n_pch[1], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[a]), fg='transparent')
)
dev.off()

## 5) Correlation between environment and host & symbiont community structure

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

# Results from given level of mutualism only (e.g., obligate: omega=1)
plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & o==0.5)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')
jit_fact = .1
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'mort_rates_RDAmean_o=0.5.pdf'), height=3, width=9)
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

xyplot(cor_a ~ factor(mra) | mrb, groups = env, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab='Relative mortality of unassociated mutualist A', ylab=expression(RDA~~R^2),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		y = as.numeric(y)
		panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
		panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
		panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(m[b]), fg='transparent'),
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



## 6) Correlation between environment and host & symbiont richness
## 7) Correlation between environment and abundance 




## Plot correlations between mutualist communities and environment
## error bars show 95th percentile from 100 different starts




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


pdf(paste0(fig_dir, 'envfilt_strength_by_type&axis_RDAmean.pdf'), height=8.5, width=11)

for(t in levels(cor_summary$topo)){
for(v in 1:2){

plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & topo==t & env==v)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

print(
xyplot(cor_a ~ factor(sigA) | sigB + envfilt, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(RDA~~R^2),main=paste('Topology =',t,', Env',v),
	panel=function(x, y, subscripts){
		x = as.numeric(x)
		panel.segments(x+jit_a[1], low95s$cor_a[subscripts], x+jit_a[1], up95s$cor_a[subscripts])
		panel.segments(x+jit_b[1], low95s$cor_b[subscripts], x+jit_b[1], up95s$cor_b[subscripts])
		panel.lines(x+jit_a[1], y, col='grey50')
		panel.lines(x+jit_a[1], means$cor_b[subscripts], col='grey50')
		panel.xyplot(x+jit_a[1], y, pch=a_pch[v], col=1)
		panel.xyplot(x+jit_b[1], means$cor_b[subscripts], pch=b_pch[v], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=c(a_pch[v], b_pch[v])), 
		text=list(c('A','B')))
)
)
}}
dev.off()

# Compare effect of topology by separating axes and partners
use_pch = c(0,1,2)
jit_fact = 0.15
jit = jit_fact*c(-1,0,1)

pdf(paste0(fig_dir, 'envfilt_strength_by_type&axis&partner_RDAmean.pdf'), height=8.5, width=11)

for(p in c('cor_a','cor_b')){
for(v in 1:2){

plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & env==v)
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

title_label = paste(ifelse(p=='cor_a', 'Host','Symbiont'), '~ Env', v)

print(
xyplot(means[,p] ~ factor(sigA) | sigB + envfilt, groups=topo, data=means, ylim = c(-.1,1),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(RDA~~R^2),main=title_label,
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit[groups[subscripts]], low95s[,p][subscripts], x+jit[groups[subscripts]], up95s[,p][subscripts])
		panel.lines(x[groups[subscripts]=='one2one']+jit[1], y[groups[subscripts]=='one2one'], col='grey50')
		panel.lines(x[groups[subscripts]=='one2many']+jit[2], y[groups[subscripts]=='one2many'], col='grey50')
		panel.lines(x[groups[subscripts]=='many2many']+jit[3], y[groups[subscripts]=='many2many'], col='grey50')
		panel.xyplot(x+jit[groups[subscripts]], y, pch=use_pch[groups[subscripts]], col=1)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), text=list(levels(means$topo)))
)
)
}}
dev.off()

# Compare effects of topology and filtering type by separating filtering strength
a_pch = c(16, 1)
b_pch = c(15, 0)
jit_fact = 0.1
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'topo_vs_envfilt_RDAmean_by_filtering_strength.pdf'), height=5, width=7.5)
for(a in 2^(-2:2)){
for(b in 2^(-2:2)){

	plot_data = subset(cor_summary, measure=='rda' & summary=='mean' & sigA==a & sigB==b)
	means = subset(plot_data, stat=='mean')
	low95s = subset(plot_data, stat=='2.5%')
	up95s = subset(plot_data, stat=='97.5%')

	title_label = bquote(sigma[A]==.(a)~~sigma[B]==.(b))
	
	print(
	xyplot(cor_a ~ topo | envfilt, groups = env, data=means, ylim = c(-.1,1),
		scales=list(alternating=1), xlab='Network Topology', ylab=expression(RDA~~R^2), main=title_label,
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
	)
}}
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
		panel.xyplot(x+jit, y, pch=use_pch[groups[subscripts]], col=1)
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
		panel.xyplot(x+jit, y, pch=use_pch[groups[subscripts]], col=1)
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
		panel.xyplot(x+jit, y, pch=use_pch[groups[subscripts]], col=1)
		panel.abline(h=100, lty=3)
	}, 
	strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
		bg='transparent', var.name=expression(sigma[b], Filtering), fg='transparent'),
	key=list(space='right', points=list(pch=use_pch), 
		text=list(levels(plot_data$topo)))
)
dev.off()

# Scaled by num. symbiont species

plot_data[plot_data$topo %in% c('many2many','one2many'),'S_b'] = plot_data[plot_data$topo %in% c('many2many','one2many'),'S_b'] / 10
plot_data[plot_data$topo == 'one2one','S_b'] = plot_data[plot_data$topo == 'one2one','S_b'] / 30
means = subset(plot_data, stat=='mean')
low95s = subset(plot_data, stat=='2.5%')
up95s = subset(plot_data, stat=='97.5%')

pdf(paste0(fig_dir, 'envfilt_strength_by_type_Sb_prop_mean.pdf'), height=8.5, width=11)
xyplot(S_b ~ factor(sigA) | sigB + envfilt, groups=topo, data=means, ylim=c(0,1.2),
	scales=list(alternating=1), xlab=expression(sigma[a]), ylab=expression(Mean~~S[B]~~(Scaled)),
	panel=function(x, y, subscripts, groups){
		x = as.numeric(x)
		panel.segments(x+jit, low95s$S_b[subscripts], x+jit, up95s$S_b[subscripts])
		panel.xyplot(x+jit, y, pch=use_pch[groups[subscripts]], col=1)
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

for(var in c('Sa','Sb','NL')){
	cor_summary[,var] = as.numeric(cor_summary[,var])
	comm_summary[,var] = as.numeric(comm_summary[,var])
}

cor_summary$a2b = cor_summary$Sa / cor_summary$Sb
comm_summary$a2b = comm_summary$Sa / comm_summary$Sb

NLratios = c(1.25,1.5,2,3)
match_NLratio = function(x){
	ds = abs(x - NLratios)
	NLratios[ds==min(ds)]
}

cor_summary$NL2a = sapply(cor_summary$NL / cor_summary$Sa, match_NLratio)
comm_summary$NL2a = sapply(comm_summary$NL / comm_summary$Sa, match_NLratio)



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


#### RUN 6: Effect of changing global abundance distribution
cor_summary$envfilt = factor(cor_summary$envfilt, levels = c('none','same','opposite','all'))
comm_summary$envfilt = factor(comm_summary$envfilt, levels = c('none','same','opposite','all'))

for(var in c('maxN','r','corrA','corrB')){
	cor_summary[,var] = as.numeric(cor_summary[,var])
	comm_summary[,var] = as.numeric(comm_summary[,var])
}


a_pch = c(16, 1)
b_pch = c(15, 0)
jit_fact = .03
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)


pdf(paste0(fig_dir, 'gsad_corr_RDAmean.pdf'), height=10, width=10)
for(f in levels(cor_summary$envfilt)){
for(n in 2^(1:5)){

	plot_data = subset(cor_summary, envfilt==f&summary=='mean'&measure=='rda'&maxN==n)

	means = subset(plot_data, stat=='mean')
	low95s = subset(plot_data, stat=='2.5%')
	up95s = subset(plot_data, stat=='97.5%')

	print(
	xyplot(cor_a ~ r | corrA + corrB, groups = env, data=means, ylim = c(0,1),
		scales=list(alternating=1), xlab='r', ylab=expression(RDA~~R^2), main=paste('Max N =',n, ', Filtering =',f),
		panel=function(x, y, subscripts, groups){
			panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
			panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
			panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
			panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
		}, 
		strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
			bg='transparent', var.name=expression(Env[A],Env[B]), fg='transparent'),
		key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
			text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
	)
	)
}}

dev.off()

combos = expand.grid(corrA=0:2, corrB=0:2)
jit_fact = .2
jit_a = -1*c(3*jit_fact/2, jit_fact/2)
jit_b = c(jit_fact/2, 3*jit_fact/2)

pdf(paste0(fig_dir, 'gsad_corr_RDAmean_compare_maxN.pdf'), height=12, width=7)
for(i in 1:nrow(combos)){

	plot_data = subset(cor_summary, corrA==combos[i,1]&corrB==combos[i,2]&summary=='mean'&measure=='rda')

	means = subset(plot_data, stat=='mean')
	low95s = subset(plot_data, stat=='2.5%')
	up95s = subset(plot_data, stat=='97.5%')

	print(
	xyplot(cor_a ~ maxN | r + envfilt, groups = env, data=means, ylim = c(0,1),
		scales=list(alternating=1), xlab='Max N', ylab=expression(RDA~~R^2), main=paste('Env A =',combos[i,1], ', Env B =',combos[i,2]),
		panel=function(x, y, subscripts, groups){
			print(groups)
			print(subscripts)
			panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
			panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
			panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
			panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
		}, 
		strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
			bg='transparent', var.name=c('r','Filtering'), fg='transparent'),
		key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
			text=list(expression(A*" ~ "*E[1],A*" ~ "*E[2],B*" ~ "*E[1],B*" ~ "*E[2])))
	)
	)

}
dev.off()

# Effect on species richness correlation
pdf(paste0(fig_dir, 'gsad_corr_Smean_compare_maxN.pdf'), height=12, width=7)
for(i in 1:nrow(combos)){

	plot_data = subset(cor_summary, corrA==combos[i,1]&corrB==combos[i,2]&summary=='mean'&measure=='S')

	means = subset(plot_data, stat=='mean')
	low95s = subset(plot_data, stat=='2.5%')
	up95s = subset(plot_data, stat=='97.5%')

	print(
	xyplot(cor_a ~ maxN | r + envfilt, groups = env, data=means, ylim = c(-1,1),
		scales=list(alternating=1), xlab='Max N', ylab='r', main=paste('Env A =',combos[i,1], ', Env B =',combos[i,2]),
		panel=function(x, y, subscripts, groups){
			panel.segments(x+jit_a[groups[subscripts]], low95s$cor_a[subscripts], x+jit_a[groups[subscripts]], up95s$cor_a[subscripts])
			panel.segments(x+jit_b[groups[subscripts]], low95s$cor_b[subscripts], x+jit_b[groups[subscripts]], up95s$cor_b[subscripts])
			panel.xyplot(x+jit_a[groups[subscripts]], y, pch=a_pch[groups[subscripts]], col=1)
			panel.xyplot(x+jit_b[groups[subscripts]], means$cor_b[subscripts], pch=b_pch[groups[subscripts]], col=1)
		}, 
		strip = strip.custom(strip.names=T, strip.levels=T, sep='=', 
			bg='transparent', var.name=c('r','Filtering'), fg='transparent'),
		key=list(space='right', points=list(pch=c(a_pch, b_pch)), 
			text=list(expression(S[A]*" ~ "*E[1],S[A]*" ~ "*E[2],S[B]*" ~ "*E[1],S[B]*" ~ "*E[2])))
	)
	)

}
dev.off()


## Plot expected number of colonists for each set of parameters


maxN_vec = 2^(1:5)
r_vec = c(0, 0.5, 0.9)
envfilt_vec = c('opposite','same','none','all')
a_corr = 0:2
b_corr = 0:2

setwd(paste0(working_dir, 'Parms/Parms_6/'))
parm_filelist = list.files(recursive=T)

# Re-order file list for plotting
runIDs = gsub('^[0-9]+-[0-9.]+/p_', '', parm_filelist)
runIDs = gsub('.txt', '', runIDs)
parm_vals = get_parms(runIDs[1])
for(i in 2:length(runIDs)){
	parm_vals = rbind(parm_vals, get_parms(runIDs[i]))
}
parm_vals$maxN = as.numeric(parm_vals$maxN)
parm_vals$r = as.numeric(parm_vals$r)
parm_vals$envfilt = factor(parm_vals$envfilt, levels=c('none','same','opposite','all'))
parm_vals$f_order = 1:nrow(parm_vals)

# Go through files and append together into a data frame

col_summary = sapply(parm_filelist, function(f){
	
	# Load this set of parameters
	source(f)
	
	# Generate random niches for each mutualist
	niche_gsad_a = make_niches_gsad(S_a, nicheparms_a, gsad_dist_a)
	niche_gsad_b = make_niches_gsad(S_b, nicheparms_b, gsad_dist_b)
	niches_a = niche_gsad_a$niches # array of S_a x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
	niches_b = niche_gsad_b$niches # array of S_b x 2 matrix of niche optima and niche breadths (mu and sigma of normal distribution)
	gsad_a = niche_gsad_a$gsad
	gsad_b = niche_gsad_b$gsad

	# Calculate expected number of colonists
	col_a = calc_colonists(niches_a, gsad_a, draw_plot=F)
	col_b = calc_colonists(niches_b, gsad_b, draw_plot=F)

	parm_vals = get_parms(runID)
	
	# Bind values to data
	c(parm_vals, list(col_a=col_a, col_b=col_b))

})
col_df = data.frame(t(col_summary))
rownames(col_df)=1:nrow(col_df)

col_df$maxN = as.numeric(unlist(col_df$maxN))
col_df$r = as.numeric(unlist(col_df$r))
col_df$corrA = as.numeric(unlist(col_df$corrA))
col_df$corrB = as.numeric(unlist(col_df$corrB))
col_df$envfilt = factor(unlist(col_df$envfilt), levels=c('none','same','opposite','all'))

col_df = col_df[with(col_df, order(maxN, corrA, corrB, r, envfilt)),]

combos = expand.grid(corrB=b_corr, corrA=a_corr, maxN=maxN_vec)
combos = subset(combos, corrB!=0|corrA!=0)

bw_col = colorRampPalette(c('white','black'))(10)
envrange=seq(-2,2,.2)
use_labs = seq(1,length(envrange), 5)

pdf(paste0(fig_dir, 'expected_num_colonists_gsad-env_correlation.pdf'), height=9, width=7)
for(i in 1:nrow(combos)){
	a = combos[i,'corrA']
	b = combos[i,'corrB']
	n = combos[i,'maxN']	

	plot_data = subset(col_df, maxN==n & corrA==a & corrB==b)
	
	layout(matrix(c(1:8, rep(9,4)), nrow=4))
	par(mar=c(.1,.1,.1,.1))
	par(oma=c(4,10,4,1))
	for(i in 1:nrow(plot_data)){
		image(plot_data$col_a[[i]], col=bw_col, axes=F)
		abline(h=0.5, v=0.5, col=2)
		box()

		if(i%%4 == 0){
			axis(1, at=(use_labs-1)/(length(envrange)-1), labels=envrange[use_labs])
			mtext('Env 1', 1, 2.5, cex=.8)
		}
		if(i <= 4){
			axis(2, at=(use_labs-1)/(length(envrange)-1), labels=envrange[use_labs], las=1)
			mtext('Env 2', 2, 2.5, cex=.8)
			mtext(plot_data[i,'envfilt'], 2, 5, las=1, adj=1)
		}
		if(i%%4 ==1){
			mtext(paste('r =',plot_data[i,'r']), 3, 1)
		}	
	}
	plot.new()
	plot.window(xlim=c(0,1), ylim=c(0,1))
	plotColorRamp(bw_col, 10, c(0,.3,.2,.7), labels=seq(0, 30, 3), 'Expected Num. Species')

	mtext(paste('MaxN =', n, ', A ~', a, ', B ~', b, ', Hosts Colonizing'), 3, 2, outer=T)

	for(i in 1:nrow(plot_data)){
		image(plot_data$col_b[[i]], col=bw_col, axes=F)
		abline(h=0.5, v=0.5, col=2)
		box()

		if(i%%4 == 0){
			axis(1, at=(use_labs-1)/(length(envrange)-1), labels=envrange[use_labs])
			mtext('Env 1', 1, 2.5, cex=.8)
		}
		if(i <= 4){
			axis(2, at=(use_labs-1)/(length(envrange)-1), labels=envrange[use_labs], las=1)
			mtext('Env 2', 2, 2.5, cex=.8)
			mtext(plot_data[i,'envfilt'], 2, 5, las=1, adj=1)
		}
		if(i%%4 ==1){
			mtext(paste('r =',plot_data[i,'r']), 3, 1)
		}	
	}
	plot.new()
	plot.window(xlim=c(0,1), ylim=c(0,1))
	plotColorRamp(bw_col, 10, c(0,.3,.2,.7), labels=seq(0, 10, 1), 'Expected Num. Species')

	mtext(paste('MaxN =', n, ', A ~', a, ', B ~', b, ', Symbionts Colonizing'), 3, 2, outer=T)
}

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

