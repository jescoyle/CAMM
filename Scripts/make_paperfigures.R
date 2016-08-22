## This script generates figures for the CAMM Manuscript


options(stringsAsFactors=F)
library(reshape2)
library(lattice) # Plots
library(latticeExtra)
library(MBESS) # confidence intervals on Cohen's d

# Main directory
working_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/'
setwd(working_dir)

# Directory for saving figures
fig_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/Figures/Manuscript/'

# Load functions
code_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/GitHub/CAMM/R/'
git_dir = 'C:/Users/jrcoyle/Documents/Research/CAMM/GitHub/CAMM/'
source(paste(code_dir,'analysis_functions.R', sep=''))
source(paste(code_dir,'simulation_functions.R', sep=''))

# Load commununity and correlation summaries from runs
comm1 = read.csv(paste0(git_dir, 'Results/comm_summary_run1.csv'))
cor1 = read.csv(paste0(git_dir, 'Results/cor_summary_run1.csv'))
comm2 = read.csv(paste0(git_dir, 'Results/comm_summary_run2.csv'))
cor2 = read.csv(paste0(git_dir, 'Results/cor_summary_run2.csv'))
comm3 = read.csv(paste0(git_dir, 'Results/comm_summary_run3.csv'))
cor3 = read.csv(paste0(git_dir, 'Results/cor_summary_run3.csv'))
comm3$o = 1 ; cor3$o = 1 # Add column indicating that this was run under obligate mutualism
comm3b = read.csv(paste0(git_dir, 'Results/comm_summary_run3b.csv'))
cor3b = read.csv(paste0(git_dir, 'Results/cor_summary_run3b.csv'))
comm3b$o = 0 ; cor3b$o = 0 # Add column indicating that this was run without mutualism



## Define symbols used for plotting throughout
a_pch = c(16, 1)
b_pch = c(15, 0)
n_pch = c(17, 2)
topo_pch = c(15,16,17); names(topo_pch) = c('one2one','one2many','many2many')

use_stats= c('S_a','S_b','Beta_a','Beta_b','N','Cor_ab')
stat_names = expression(S[host],S[sym],beta[host],beta[sym],N,S[host]%~%S[sym]); names(stat_names)=use_stats
stat_pch = c(0,1,2,4,5,6); names(stat_pch) = use_stats

use_cor = c('cor_a_1','cor_a_2','cor_b_1','cor_b_2')
cor_names = expression(Host%~%E[1], Host%~%E[2], Sym%~%E[1], Sym%~%E[2]); names(cor_names)=use_cor
cor_S_names = expression(S[host]%~%E[1], S[host]%~%E[2], S[sym]%~%E[1], S[sym]%~%E[2])
cor_pch = c(a_pch, b_pch); names(cor_pch) = use_cor

envfilt = c('none','same','opposite','all')
env_names = c('None','Same Gradient','Opposite Gradients','Both Gradients'); names(env_names) = c('none','same','opposite','all')
env_names_stacked = c('None','Same\nGradient','Opposite\nGradients','Both\nGradients'); names(env_names_stacked) = c('none','same','opposite','all')
topos = c('one2one','one2many','many2many')
topo_names = c('Both Specialist', 'Host Specialist', 'Both Non-Specialist'); names(topo_names) = c('one2one','one2many','many2many')

##################################################################
## Functions

# A function that calculates standardized mean difference (Cohen's d)
# m : a vector of length 2 with the sample means
# s : a vector of length 2 with the sample standard deviations
# n : a vector of length 2 with the sample sizes
# use.s : character indicating whether the standard deviation used should be 'pooled', 's1', or 's2'

calc_d = function(m, s, n, use.s='pooled'){

	if(use.s=='pooled'){
		s = sqrt(((n[1]-1)*s[1]^2 + (n[2]-1)*s[2]^2)/(sum(n)+2))	
	}
	if(use.s=='s1'){
		s = s[1]
	}
	if(use.s=='s2'){
		s = s[2]
	}
	
	(m[2]- m[1])/s	
}


# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels  : vector of labels for bar
# title   : title to print above bar
# mycex   : size of label and title text
# ndig    : number of digits to round numeric labels
plotColorBoxes = function(cols, barends, labels, title=NA, mycex=1.5, ndig=0){
	n = length(cols)
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	lend.old = par('lend')
	par(lend=1)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=cols[i], border=NA)
	}

	if(is.numeric(labels)){
		labels = format(round(labels, ndig), nsmall=ndig, trim=F)
	}


	Yposition = barends[2] + dy*(0:n)

	text(barends[3]+dX*0.5, Yposition, labels, pos=4, cex=mycex)
	segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	

	if(!is.na(title)){
		
		## Determine how many characters away to place title
		digits = max(nchar(labels)) # Maximum number of digits in a label
		largest = labels[which(nchar(labels)==digits)] # Which labels are longest
		
		small.chars = grep('[-.]', largest) # Does the largest label have a small character?
			if(length(small.chars)==length(largest)) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
	par(lend=lend.old)
}



##################################################################
## Experiment 1: Degree of Mutualistic Dependence

## Figure 1
## To show no effect of weak facultative mutualism: Effects size plot of 
## (A) all community stats and (B) community correlations from omega=0 baseline (from RUN1 with mort rates at 10) regressed on omega.

omegas = seq(0,1,.1)

# A) Community statistics
plot_data = subset(comm1, mra==10 & mrb==10 & summary=='mean')
plot_data = melt(plot_data, id.vars=c('o','stat'), measure.vars=6:14)
plot_data = acast(plot_data, o~stat~variable)

# Calculate standardized mean difference (Cohen's d)
# Use both pooled and s1 for variance
eff_data = sapply(use_stats, function(y){
	as.numeric(sapply(omegas, function(o){
		calc_d(plot_data[c('0',o),'mean',y], sqrt(plot_data[c('0',o),'var',y]), c(100,100), use.s='pooled')
	}))
})

# Calculate CI on Cohen's d (standardized mean difference), assumes equal variance
CIs = apply(eff_data, c(1,2), function(x){
	ci = ci.smd(smd=x, n.1=100, n.2=100, conf.level=0.95)
	as.numeric(ci)
})

#svg(paste0(fig_dir,'Effect of mutualism strength on community statistics pooled.svg'), height=4.5, width=5)
par(mar=c(5,4,1,4))
par(lend=1)
plot(omegas, rep(0,11), type='n', ylim=range(CIs), axes=F,
	xlab=expression(Strength~~of~~mutualism~~(omega)), ylab='Effects Size (d)')
abline(h=0, col='grey50', lty=2, lend=1)
for(stat in use_stats){
	polygon(c(omegas,rev(omegas)), c(CIs[1,,stat],CIs[3,11:1,stat]), col='#00000020', border=NA)	
	points(omegas, eff_data[,stat], pch=stat_pch[stat])
	lines(omegas, eff_data[,stat])
}
par(xpd=T)
text(rep(1, length(use_stats)), eff_data[11,], stat_names, pos=4)
par(xpd=F)
#legend('bottomleft', stat_names, pch=use_pch, bty='n')
axis(1)
axis(2, las=1)
abline(h=par('usr')[3], lwd=3)
abline(v=par('usr')[1], lwd=3)
#dev.off()

# B) Correlations
plot_data = subset(cor1, mra==10 & mrb==10 & summary=='mean')
plot_data = melt(plot_data, id.vars=c('o','stat','env','measure'), measure.vars=c('cor_a','cor_b'))
plot_data = acast(plot_data, o ~ stat ~ measure ~ variable + env)

# Calculate standardized mean difference (Cohen's d)
# Use both pooled and s1 for variance
eff_data = apply(plot_data, c(3,4), function(y){
	as.numeric(sapply(omegas, function(o){
		calc_d(y[c('0',o),'mean'], sqrt(y[c('0',o),'var']), c(100,100), use.s='s1')
	}))
})

# Calculate CI on Cohen's d (standardized mean difference), assumes equal variance
CIs = apply(eff_data, c(1,2,3), function(x){
	ci = ci.smd(smd=x, n.1=100, n.2=100, conf.level=0.95)
	as.numeric(ci)
})


# Community composition
#svg(paste0(fig_dir,'Effect of mutualism strength on community composition correlation s1.svg'), height=4.5, width=5)
par(mar=c(5,4,1,4))
par(lend=1)
plot(omegas, rep(0,11), type='n', ylim=range(CIs[,,'rda',]), axes=F,
	xlab=expression(Strength~~of~~mutualism~~(omega)), ylab='Std. mean difference (d)')
abline(h=0, col='grey50', lty=2, lend=1)
for(stat in use_cor){
	polygon(c(omegas,rev(omegas)), c(CIs[1,,'rda',stat],CIs[3,11:1,'rda',stat]), col='#00000020', border=NA)	
	points(omegas, eff_data[,'rda',stat], pch=cor_pch[stat])
	lines(omegas, eff_data[,'rda',stat])
}
legend('topleft', cor_names, pch=cor_pch, bty='n', title=expression(RDA~~R^2), ncol=2)
axis(1)
axis(2, las=1)
abline(h=par('usr')[3], lwd=3)
abline(v=par('usr')[1], lwd=3)
#dev.off()

# Species richness
#svg(paste0(fig_dir,'Effect of mutualism strength on species richness correlation s1.svg'), height=4.5, width=5)
par(mar=c(5,4,1,4))
par(lend=1)
plot(omegas, rep(0,11), type='n', ylim=range(CIs[,,c('rda','S'),]), axes=F,
	xlab=expression(Strength~~of~~mutualism~~(omega)), ylab='Std. mean difference (d)')
abline(h=0, col='grey50', lty=2, lend=1)
for(stat in use_cor){
	polygon(c(omegas,rev(omegas)), c(CIs[1,,'S',stat],CIs[3,11:1,'S',stat]), col='#00000020', border=NA)	
	points(omegas, eff_data[,'S',stat], pch=cor_pch[stat])
	lines(omegas, eff_data[,'S',stat])
}
legend('topleft', cor_S_names, pch=cor_pch, bty='n', title=expression(Pearson~~italic(r)), ncol=2)
axis(1)
axis(2, las=1)
abline(h=par('usr')[3], lwd=3)
abline(v=par('usr')[1], lwd=3)
#dev.off()


## Figure 2
## To show interaction with topology and type of filtering: Effects size (from o=0) of 
## (A) (some) community stats and (B) community correlations under different topology and filtering (from RUN 2).

omegas = c(0,0.5,0.9,1)

## A) Community statistics
plot_data = subset(comm2, summary=='mean')
plot_data = melt(plot_data, id.vars=c('o','topo','envfilt','stat'), measure.vars=6:14)
plot_data = acast(plot_data, o~topo~envfilt~stat~variable)

# Calculate standardized mean difference (Cohen's d)
# Use both pooled and s1 for variance
eff_data = apply(plot_data[,,,,use_stats], c(2,3,5), function(y){
	as.numeric(sapply(rownames(y), function(o){
		calc_d(y[c('0',o),'mean'], sqrt(y[c('0',o),'var']), c(100,100), use.s='pooled')
	}))
})

# Calculate CI on Cohen's d (standardized mean difference), assumes equal variance
CIs = apply(eff_data, c(1,2,3,4), function(x){
	if(is.na(x)){
		rep(NA,3)
	} else {
		as.numeric(ci.smd(smd=x, n.1=100, n.2=100, conf.level=0.95))
	}
})


# Cor a ~ b is NA for topo = many2many and envfilt = none because S_b = 10 in all communities for several chains
# May want to re-write summary code so that NAs don't propagate, but will wait to see how pervasive this problem is.

#svg(paste0(fig_dir, 'Effect of mutualism strength on community stats across topo and envfilt (pooled).svg'))
par(mar=c(.5,.5,.5,.5))
par(oma=c(4,4,3,0))
layout(matrix(c(1:12, rep(13,4)), nrow=4))
par(lend=1)

for(n in topos){
for(f in envfilt){

plot(omegas, rep(0,4), type='n', ylim=range(CIs, na.rm=T), axes=F,
	xlab='', ylab='')
abline(h=0, col='grey50', lty=2, lend=1)
for(stat in use_stats){
	not_missing = which(!is.na(eff_data[,n,f,stat]))
	use_x = omegas[not_missing]
	polygon(c(use_x,rev(use_x)), c(CIs[1,not_missing,n,f,stat],CIs[3,rev(not_missing),n,f,stat]), col='#00000020', border=NA)	
	points(use_x, eff_data[not_missing,n,f,stat], pch=stat_pch[stat])
	lines(use_x, eff_data[not_missing,n,f,stat])
}
# Add labels
if(f == 'all') axis(1)
if(n == 'one2one') axis(2, las=1)
if(f == 'none') mtext(topo_names[n], 3, 1, cex=.8)
if(n == 'many2many'){
	par(xpd=NA)
	text(1, 1 , env_names[f], srt=-90, cex=1.3, adj=c(.5,-2))
	par(xpd=F)
}
	
# Add axis lines
abline(h=par('usr')[3], lwd=3)
abline(v=par('usr')[1], lwd=3)


}} # closes for loops

mtext(expression(Strength~~of~~mutualism~~(omega)), 1, 2.5, cex=.8, outer=T, adj=0.35)
mtext('Std. mean difference', 2, 2.5, cex=.8, outer=T)

plot.new()
legend('center', stat_names, pch=use_pch, bty='n', cex=1.2)

#dev.off()


## Plot only for obligate mutualism

jit = c(-.1,0,.1); names(jit) = topos
x = 1:length(use_stats)

#svg(paste0(fig_dir,'Effects of obligate mutualism on community stats across envfilt and topo (pooled).svg'), height=9, width=3.5)
par(mfrow=c(4,1))
par(mar=c(3,4,2,.5))
par(lend=1)
for(f in envfilt){
	
	plot(x, rep(0, length(use_stats)), type='n', xlim=c(0.5,length(use_stats)+.5),
		ylim=range(CIs[,4,,,], na.rm=T),axes=F, ylab='Std. mean difference', xlab='')
	
	# Add axis lines
	axis(1, at=x, labels=stat_names, tick=F, padj=0)
	axis(2, las=1)
	#abline(v=x+.5, lty=1, col='grey80')
	abline(h=seq(-4,2,2), lty=2, col='grey')
	abline(h=par('usr')[3], lwd=3)
	abline(v=par('usr')[1], lwd=3)

	# Add points
	for(n in topos){
		xvals = jit[n] + x
		arrows(xvals, CIs[1,4,n,f,], xvals, CIs[3,4,n,f,], angle=90, length=.03, code=3)
		points(xvals, eff_data[4,n,f,use_stats], pch=topo_pch[n])
	}

	# Add plot label
	plot_lab = paste0(letters[which(envfilt==f)], ') ', env_names[f])
	mtext(plot_lab, 3, at=-1, line=0.5, adj=0, cex=.8)

	# Add legend to first plot
	if(f == 'none') legend('topright', topo_names[topos], pch=topo_pch[topos], border='grey', bg='white')

}
#dev.off()


## B) Community correlations
plot_data = subset(cor2, summary=='mean')
plot_data = melt(plot_data, id.vars=c('o','topo','envfilt','stat','env','measure'), measure.vars=c('cor_a','cor_b'))
plot_data = acast(plot_data, o ~ topo~envfilt ~ stat ~ measure ~ variable + env)

# Calculate standardized mean difference (Cohen's d)
# Use both pooled and s1 for variance
eff_data = apply(plot_data, c(2,3,5,6), function(y){
	as.numeric(sapply(omegas, function(o){
		calc_d(y[c('0',o),'mean'], sqrt(y[c('0',o),'var']), c(100,100), use.s='pooled')
	}))
})

# Calculate CI on Cohen's d (standardized mean difference), assumes equal variance
CIs = apply(eff_data, 1:5, function(x){
	if(is.na(x)){
		rep(NA,3)
	} else {
		as.numeric(ci.smd(smd=x, n.1=100, n.2=100, conf.level=0.95))
	}
})

svg(paste0(fig_dir, 'Effect of mutualism strength on community correlations across topo and envfilt (pooled).svg'))
par(mar=c(.5,.5,.5,.5))
par(oma=c(4,4,3,0))
layout(matrix(c(1:12, rep(13,4)), nrow=4))
par(lend=1)

for(n in topos){
for(f in envfilt){

plot(omegas, rep(0,4), type='n', ylim=range(CIs[,,,,'rda',], na.rm=T), axes=F,
	xlab='', ylab='')
abline(h=0, col='grey50', lty=2, lend=1)
for(stat in use_cor){
	not_missing = which(!is.na(eff_data[,n,f,'rda',stat]))
	use_x = omegas[not_missing]
	polygon(c(use_x,rev(use_x)), c(CIs[1,not_missing,n,f,'rda',stat],CIs[3,rev(not_missing),n,f,'rda',stat]), col='#00000020', border=NA)	
	points(use_x, eff_data[not_missing,n,f,'rda',stat], pch=cor_pch[stat])
	lines(use_x, eff_data[not_missing,n,f,'rda',stat])
}
# Add labels
if(f == 'all') axis(1)
if(n == 'one2one') axis(2, las=1)
if(f == 'none') mtext(topo_names[n], 3, 1, cex=.8)
if(n == 'many2many'){
	par(xpd=NA)
	text(1, 1 , env_names[f], srt=-90, cex=1.3, adj=c(.5,-2))
	par(xpd=F)
}
	
# Add axis lines
abline(h=par('usr')[3], lwd=3)
abline(v=par('usr')[1], lwd=3)

}} # closes for loops

mtext(expression(Strength~~of~~mutualism~~(omega)), 1, 2.5, cex=.8, outer=T, adj=0.35)
mtext('Std. mean difference', 2, 2.5, cex=.8, outer=T)

plot.new()
legend('center', cor_names, pch=cor_pch, bty='n', cex=1.2, title=expression(RDA~~R^2))

dev.off()

# Plot only for obligate mutualism
jit = c(-.1,0,.1); names(jit) = topos
x = 1:length(use_cor)

svg(paste0(fig_dir,'Effects of obligate mutualism on community correlations across envfilt and topo (pooled).svg'), height=9, width=3.5)
par(mfrow=c(4,1))
par(mar=c(3,4,2,.5))
par(lend=1)
for(f in envfilt){
	
	plot(x, rep(0, length(use_cor)), type='n', xlim=c(0.5,length(use_cor)+.5),
		ylim=range(CIs[,4,,,'rda',], na.rm=T),axes=F, ylab='Std. mean difference', xlab='')
	
	# Add axis lines
	axis(1, at=x, labels=cor_names, tick=F, padj=0)
	axis(2, las=1)
	#abline(v=x+.5, lty=1, col='grey80')
	abline(h=seq(-2,4,2), lty=2, col='grey')
	abline(h=par('usr')[3], lwd=3)
	abline(v=par('usr')[1], lwd=3)

	# Add points
	for(n in topos){
		xvals = jit[n] + x
		arrows(xvals, CIs[1,4,n,f,'rda',use_cor], xvals, CIs[3,4,n,f,'rda',use_cor], angle=90, length=.03, code=3)
		points(xvals, eff_data[4,n,f,'rda',use_cor], pch=topo_pch[n])
	}

	# Add plot label
	plot_lab = paste0(letters[which(envfilt==f)], ') ', env_names[f])
	mtext(plot_lab, 3, at=-.5, line=0.5, adj=0, cex=.8)

	# Add legend to first plot
	if(f == 'none') legend('topright', topo_names[topos], pch=topo_pch[topos], border='grey', bg='white')

}
dev.off()


##################################################################
## Experiment 2: Strength and Type of Environmental Constraint

# Combine results from run 3 with and without mutualism
comm3 = rbind(comm3, comm3b)
cor3 = rbind(cor3, cor3b)

## TWO WAYS TO DISPLAY: AS EFFECTS SIZE OF MUTUALISM OR AS TWO SEPARATE SIMULATIONS W/ AND W/O MUTUALISM


# A) Community statistics
plot_data = subset(comm3, summary=='mean')
plot_data = melt(plot_data, id.vars=c('o','sigA','sigB','topo','envfilt','stat'), measure.vars=7:15)
plot_data = acast(plot_data, o~topo~envfilt~sigA~sigB~stat~variable)

# Calculate standardized mean difference (Cohen's d)
# Use both pooled and s1 for variance
# Effects for S_tot are NA when variance is 0 b/c all species always present
eff_data = apply(plot_data[,,,,,,use_stats], c(2,3,4,5,7), function(X){
	calc_d(X[,'mean'], sqrt(X[,'var']), c(100,100), use.s='pooled')
})


# Calculate CI on Cohen's d (standardized mean difference), assumes equal variance
CIs = apply(eff_data, 1:5, function(x){
	if(is.na(x)){
		rep(NA,3)
	} else {
		as.numeric(ci.smd(smd=x, n.1=100, n.2=100, conf.level=0.95))
	}
})


# Display effects as heat map

use_col = colorRampPalette(c('darkblue','white','darkred'))(10)

# Each network topology on separate pages
for(n in topos){

data_range = range(eff_data[n,,,,], na.rm=T)
data_range[1] = min(floor(data_range[1]),-8)
data_range[2] = max(ceiling(data_range[2]),8)
use_breaks = c(data_range[1],-8,-4,-2,-1,0,1,2,4,8, data_range[2])

#pdf(paste0(fig_dir, 'Effects of filtering strength on community stats ', n, '.pdf'), height=9, width=7)
	par(lend=1)
	layout(matrix(c(1:18, rep(19,6)), nrow=6))
	par(mar=c(.25,.25,.25,.25))
	par(oma=c(4.5,9,2,0))

	for(f in envfilt[2:4]){
	for(y in use_stats){
		use_data = eff_data[n,f,,,y]
		image(1:5,1:5, use_data, breaks=use_breaks, col=use_col, axes=F, xlab='', ylab='')

		missing = which(is.na(use_data), arr.ind=T)
		abline(h=c(1.5, 2.5, 3.5, 4.5), col='white')
		abline(v=c(1.5, 2.5, 3.5, 4.5), col='white')
		points(missing[,1], missing[,2],pch=4, col='grey', cex=1.5)

		if(y==use_stats[6]){
			axis(1, at=1:5, labels=rownames(use_data), las=2, tick=F, line=0)
			mtext(expression(sigma[host]^2), 1, 3.5, cex=.8)
		}
		
		if(y==use_stats[1]){
			mtext(env_names[f], 3, 0.5)
		}

		if(f==envfilt[2]){

			axis(2, at=1:5, labels=colnames(use_data), las=1, tick=F, line=0)
			
			par(xpd=NA)
			text(-3.7, 5, stat_names[y], pos=4, cex=1.5)
			par(xpd=F)

			mtext(expression(sigma[sym]^2), 2, 3, cex=.8)
		}
	
	}} # closes envfilt and use stats loops

	plot.new()
	plotColorBoxes(use_col, c(.1,.2,.2,.7), use_breaks, title='Std. mean difference')

#dev.off()
} # closes topo loop


# Display only host statistics since symbionts not expected to differ from no mutualism situation
host_stats = c('S_a','Beta_a','N','Cor_ab')

for(n in topos){

data_range = range(eff_data[n,,,,host_stats], na.rm=T)
data_range[1] = min(floor(data_range[1]),-8)
data_range[2] = max(ceiling(data_range[2]),8)
use_breaks = c(data_range[1],-8,-4,-2,-1,0,1,2,4,8, data_range[2])

#pdf(paste0(fig_dir, 'Effects of filtering strength on community stats ', n, ' host.pdf'), height=6, width=7)
	par(lend=1)
	layout(matrix(c(1:12, rep(13,4)), nrow=4))
	par(mar=c(.25,.25,.25,.25))
	par(oma=c(4.5,9,2,0))

	for(f in envfilt[2:4]){
	for(y in host_stats){
		use_data = eff_data[n,f,,,y]
		image(1:5,1:5, use_data, breaks=use_breaks, col=use_col, axes=F, xlab='', ylab='')

		missing = which(is.na(use_data), arr.ind=T)
		abline(h=c(1.5, 2.5, 3.5, 4.5), col='white')
		abline(v=c(1.5, 2.5, 3.5, 4.5), col='white')
		points(missing[,1], missing[,2],pch=4, col='grey', cex=1.5)

		if(y==host_stats[4]){
			axis(1, at=1:5, labels=rownames(use_data), las=2, tick=F, line=0)
			mtext(expression(sigma[host]^2), 1, 3.5, cex=.8)
		}
		
		if(y==host_stats[1]){
			mtext(env_names[f], 3, 0.5)
		}

		if(f==envfilt[2]){

			axis(2, at=1:5, labels=colnames(use_data), las=1, tick=F, line=0)
			
			par(xpd=NA)
			text(-3.7, 5, stat_names[y], pos=4, cex=1.5)
			par(xpd=F)

			mtext(expression(sigma[sym]^2), 2, 3, cex=.8)
		}
	
	}} # closes envfilt and use stats loops

	plot.new()
	plotColorBoxes(use_col, c(.1,.2,.2,.7), use_breaks, title='Std. mean difference')

#dev.off()
} # closes topo loop

# B) Environmental correlations
plot_data = subset(cor3, summary=='mean'&measure=='rda')
plot_data = melt(plot_data, id.vars=c('o','sigA','sigB','topo','envfilt','stat','env'), measure.vars=c('cor_a','cor_b'))
plot_data = acast(plot_data, o~topo~envfilt~sigA~sigB~stat~variable+env)

# Calculate standardized mean difference (Cohen's d)
# Use both pooled and s1 for variance
# Effects for S_tot are NA when variance is 0 b/c all species always present
eff_data = apply(plot_data, c(2,3,4,5,7), function(X){
	calc_d(X[,'mean'], sqrt(X[,'var']), c(100,100), use.s='pooled')
})

# Calculate CI on Cohen's d (standardized mean difference), assumes equal variance
CIs = apply(eff_data, 1:5, function(x){
	if(is.na(x)){
		rep(NA,3)
	} else {
		as.numeric(ci.smd(smd=x, n.1=100, n.2=100, conf.level=0.95))
	}
})

# Display effects as heat map

use_col = colorRampPalette(c('darkblue','white','darkred'))(10)

# Each network topology on separate pages
for(n in topos){

#data_range = range(eff_data[n,,,,], na.rm=T)
#data_range[1] = min(floor(data_range[1]),-8)
#data_range[2] = max(ceiling(data_range[2]),8)
#use_breaks = c(data_range[1],-8,-4,-2,-1,0,1,2,4,8, data_range[2])
use_breaks = c(-8,-6,-4,-2,-1,0,1,2,4,6,8)

#pdf(paste0(fig_dir, 'Effects of filtering strength on rda R2 ', n, '.pdf'), height=6, width=7)
	par(lend=1)
	layout(matrix(c(1:12, rep(13,4)), nrow=4))
	par(mar=c(.25,.25,.25,.25))
	par(oma=c(4.5,9,2,0))

	for(f in envfilt[2:4]){
	for(y in use_cor){
		use_data = eff_data[n,f,,,y]
		image(1:5,1:5, use_data, breaks=use_breaks, col=use_col, axes=F, xlab='', ylab='')

		missing = which(is.na(use_data), arr.ind=T)
		abline(h=c(1.5, 2.5, 3.5, 4.5), col='white')
		abline(v=c(1.5, 2.5, 3.5, 4.5), col='white')
		points(missing[,1], missing[,2],pch=4, col='grey', cex=1.5)

		if(y==use_cor[4]){
			axis(1, at=1:5, labels=rownames(use_data), las=2, tick=F, line=0)
			mtext(expression(sigma[host]^2), 1, 3.5, cex=.8)
		}
		
		if(y==use_cor[1]){
			mtext(env_names[f], 3, 0.5)
		}

		if(f==envfilt[2]){

			axis(2, at=1:5, labels=colnames(use_data), las=1, tick=F, line=0)
			
			par(xpd=NA)
			text(-3.7, 5, cor_names[y], pos=4, cex=1.5)
			par(xpd=F)

			mtext(expression(sigma[sym]^2), 2, 3, cex=.8)
		}
	
	}} # closes envfilt and use stats loops

	plot.new()
	plotColorBoxes(use_col, c(.1,.2,.2,.7), use_breaks, title='Std. mean difference')

#dev.off()
} # closes topo loop

# Display for hosts only, since symbionts aren't really expected to diverge from no-mutualism situation
host_cor = c('cor_a_1','cor_a_2')
use_breaks = c(-8,-6,-4,-2,-1,0,1,2,4,6,8)

for(n in topos){
#pdf(paste0(fig_dir, 'Effects of filtering strength on rda R2 ', n, ' host.pdf'), height=4, width=7)
	par(lend=1)
	layout(matrix(c(1:6, rep(7,2)), nrow=2))
	par(mar=c(.25,.25,.25,.25))
	par(oma=c(4.5,9,2,0))

	for(f in envfilt[2:4]){
	for(y in host_cor){
		use_data = eff_data[n,f,,,y]
		image(1:5,1:5, use_data, breaks=use_breaks, col=use_col, axes=F, xlab='', ylab='')

		missing = which(is.na(use_data), arr.ind=T)
		abline(h=c(1.5, 2.5, 3.5, 4.5), col='white')
		abline(v=c(1.5, 2.5, 3.5, 4.5), col='white')
		points(missing[,1], missing[,2],pch=4, col='grey', cex=1.5)

		if(y==host_cor[2]){
			axis(1, at=1:5, labels=rownames(use_data), las=2, tick=F, line=0)
			mtext(expression(sigma[host]^2), 1, 3.5, cex=.8)
		}
		
		if(y==host_cor[1]){
			mtext(env_names[f], 3, 0.5)
		}

		if(f==envfilt[2]){

			axis(2, at=1:5, labels=colnames(use_data), las=1, tick=F, line=0)
			
			par(xpd=NA)
			text(-3.7, 5, cor_names[y], pos=4, cex=1.5)
			par(xpd=F)

			mtext(expression(sigma[sym]^2), 2, 3, cex=.8)
		}
	
	}} # closes envfilt and use stats loops

	plot.new()
	plotColorBoxes(use_col, c(.1,.2,.2,.8), use_breaks, title='Std. mean difference')

#dev.off()
} # closes topo loop


# C) Display actual values
comm_data = subset(comm3, summary=='mean')
comm_data = melt(comm_data, id.vars=c('o','sigA','sigB','topo','envfilt','stat'), measure.vars=7:15)
cor_data = subset(cor3, summary=='mean'&measure=='rda')
cor_data = melt(cor_data, id.vars=c('o','sigA','sigB','topo','envfilt','stat','env'), measure.vars=c('cor_a','cor_b'))
cor_data = dcast(cor_data, o+sigA+sigB+topo+envfilt+stat~variable+env)
cor_data = melt(cor_data, id.vars=c('o','sigA','sigB','topo','envfilt','stat'), measure.vars=7:10)
plot_data = rbind(comm_data, cor_data)
plot_data = acast(plot_data, o~topo~envfilt~sigA~sigB~stat~variable)

plot_vars = c(use_stats, use_cor)
plot_var_names = c(stat_names, cor_names)

pch0 = topo_pch; col0 = 'grey'
pch1 = 0:2; col1 = 'black'

xvals = sapply(1:5, function(x) x+c(-.25,0,.25))

#pdf(paste0(fig_dir, 'Effects of env filtering strength and type on community stats and env correlations.pdf'), height=9, width=8)

# Each variable on a separate page
for(y in plot_vars){

	par(lend=1)
	layout(matrix(c(5:1,10:6,15:11, rep(16,5)), nrow=5))
	par(mar=c(.5,.5,.5,.5))
	par(oma=c(4.5,9,4,0))

	data_range = range(plot_data[,,envfilt[2:4],,,c('2.5%','97.5%'),y], na.rm=T)
	data_range = data_range + 0.01*abs(diff(data_range))*c(-1,1)

	# Type of filtering in columns
	for(f in envfilt[2:4]){
	
	# Filtering on hosts in rows (sigA)
	for(s in dimnames(plot_data)[[4]]){
		
		# Set up plot
		frame()
		plot.window(xlim=c(0.5,5.5), ylim=data_range, xlab='', ylab='')
		
		
		if(f==envfilt[2]){
			yax = axis(2, las=1)
			if(y %in% use_stats) mtext(stat_names[y], 2, 3, cex=.8)
			if(y %in% use_cor) mtext(expression(RDA~~R^2), 2, 3)

			par(xpd=NA)
			text(-4, par('usr')[4], bquote(sigma[host]^2==.(s)), pos=4, cex=1.5)
			par(xpd=F)
		}
	
		if(s=='0.25'){
			axis(1, at=1:5, labels=dimnames(plot_data)[[5]])
			mtext(expression(sigma[sym]^2), 1, 3, cex=.8)
		}

		if(s=='4'){
			mtext(env_names[f], 3, .5)
		}

		# Add axes and labels
		abline(h=par('usr')[3], lwd=2)
		abline(h=yax, lwd=1, col='grey90')
		abline(v=par('usr')[1], lwd=2)

		# Estimates under no mutualism
		points(xvals, plot_data['0',topos,f,s,,'mean',y], pch=pch0, col=col0)
		segments(xvals, plot_data['0',topos,f,s,,'2.5%',y], xvals, plot_data['0',topos,f,s,,'97.5%',y], col=col0, lwd=2.5)
		
		# Estimates under obligate mutualism
		points(xvals, plot_data['1',topos,f,s,,'mean',y], pch=pch1, col=col1)
		segments(xvals, plot_data['1',topos,f,s,,'2.5%',y], xvals, plot_data['1',topos,f,s,,'97.5%',y], col=col1, lwd=1)
		
		
	}} # closes loops through envfilt and sigA

	# Add legend
	frame()

	points(rep(.1,3), c(.55,.5,.45), pch=pch0, col=col0)
	points(rep(.2,3), c(.55,.5,.45), pch=pch1, col=col1)
	text(rep(.25,3), c(.55,.5,.45), topo_names, pos=4)
	text(c(.1,.2), rep(.57,2), c('No Mutualism','Obligate Mutualism'), srt=90, adj=-.1)

	# Add title
	text(.5, .8, plot_var_names[y], cex=2)
}

#dev.off()


##########################################################
### Experiment 3: Effect of network topology

comm_data = subset(comm3, summary=='mean'&sigA==0.5&sigB==0.5)
comm_data = melt(comm_data, id.vars=c('o','topo','envfilt','stat'), measure.vars=7:15)
cor_data = subset(cor3, summary=='mean'& measure=='rda'&sigA==0.5&sigB==0.5)
cor_data = melt(cor_data, id.vars=c('o','topo','envfilt','stat','env'), measure.vars=c('cor_a','cor_b'))
cor_data = dcast(cor_data, o+topo+envfilt+stat~variable+env)
cor_data = melt(cor_data, id.vars=c('o','topo','envfilt','stat'), measure.vars=5:8)
plot_data = rbind(comm_data, cor_data)
plot_data = acast(plot_data, o~topo~envfilt~stat~variable)

# Four panels (envfilt)
# RDA R2 for host and sym vs. network topology (use sigmas = 0.5) 
# Show both obligate and no mutualism

pch0 = topo_pch; col0 = 'grey'
pch1 = 0:2; col1 = 'black'

plot_vars = c(use_stats, use_cor)
plot_var_names = c(stat_names, cor_names)
plot_var_units = expression(Num.~~species, Num.~~species, ' ', ' ', Num.~~individuals, italic(r), RDA~~R^2, RDA~~R^2, RDA~~R^2, RDA~~R^2)
names(plot_var_units) = plot_vars
n_y = length(plot_vars)
xvals = sapply(1:length(envfilt), function(x) x+c(-.15,0,.15))

pdf(paste0(fig_dir, 'Effects of network topology and type of filtering on community stats and env correlations.pdf'), height=9.5, width=5.7)

par(lend=1)
layout(matrix(c(1:n_y, rep(n_y+1,n_y)),nrow=n_y), widths=c(.6,.4))
par(mar=c(.65,.5,.65,.5))
par(oma=c(4,8,0,0))

for(y in plot_vars){
	data_range = range(plot_data[c('0','1'),,,c('2.5%','97.5%'),y], na.rm=T)	
	data_range = abs(diff(data_range))*0.05*c(-1,1) + data_range

	frame()
	plot.window(ylim=data_range, xlim=c(0.5, length(envfilt)+.5))
	
	# Add axes and labels
	yax = axis(2, las=1)
	abline(h=yax, col='grey90')
	abline(h=par('usr')[3], lwd=2)
	abline(v=par('usr')[1], lwd=2)
	if(y == plot_vars[n_y]) axis(1, at=xvals[2,], labels=env_names_stacked[envfilt], las=1, padj=1)
	#mtext(plot_var_units[y], 2, 3, cex=.8)
	
	this_letter = paste0(letters[which(plot_vars==y)],')')
	this_lab = plot_var_names[y][[1]]

	par(xpd=NA)
	text(-1.5, par('usr')[4], bquote(.(this_letter)~~.(this_lab)), pos=4)
	par(xpd=F)

	# Add no-mutualism points
	points(xvals, plot_data['0',topos,envfilt,'mean',y], pch=pch0, col=col0)
	segments(xvals, plot_data['0',topos,envfilt,'2.5%',y], xvals, plot_data['0',topos,envfilt,'97.5%',y], lwd=2.5, col=col0)
	
	# Add  obligate mutualism points
	points(xvals, plot_data['1',topos,envfilt,'mean',y], pch=pch1, col=col1)
	segments(xvals, plot_data['1',topos,envfilt,'2.5%',y], xvals, plot_data['1',topos,envfilt,'97.5%',y], lwd=1, col=col1)
	
}

# Add legend
frame()
points(rep(.05,3), c(.52,.5,.48), pch=pch0, col=col0)
points(rep(.15,3), c(.52,.5,.48), pch=pch1, col=col1)
text(rep(.25,3), c(.52,.5,.48), topo_names, adj=0)
text(c(.05,.15), rep(.54,2), c('No Mutualism','Obligate Mutualism'), srt=90, adj=0)

dev.off()


######################################################3
### Comparison
















