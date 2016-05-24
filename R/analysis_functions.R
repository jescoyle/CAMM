## This script contains functions for visualizing, summarizing and analyzing simulation results


### Functions for compiling multiple simulation runs ###

# A function that extracts parameter values from a runID
#	runID = a string in the filename that gives the parameters
#	assoc_str = character pairing a parameter name with its value
#	sep_str = character separating different parameter-value pairs
get_parms = function(runID, assoc_str='-', sep_str='_'){
	parm_pairs = sapply(strsplit(runID, sep_str), function(x) strsplit(x, assoc_str))
	
	vals = sapply(parm_pairs, function(x) x[2])
	vals = as.data.frame(t(vals))
	names(vals) = sapply(parm_pairs, function(x) x[1])

	vals
}


 

### Functions for visualizing simulations ###

# A function that plots the probability density functions of niches used in a simulation
#	niches = an S x 2 x 2 array giving mu and sigma parameters of the gaussian distribution for 2 environmental variables and S species
# 	grad = a 2 x 2 matrix giving the length of each environmental gradient
#	add_env = an optional matrix of environmental values at sites
plot_niches = function(niches, grad, add_env=NULL){

	par(mfrow=c(1,2))
	
	# Find number of species
	N_S = dim(niches)[1]

	# Choose colors for each species
	cols = colorRampPalette(c('#90ABFC','#000000'))(N_S)
	cols = cols[rank(niches[,'mu',1])]
	
	# Make plot for env var 1
	plot(c(0,0), xlim=grad[,1], ylim=c(0,1), xlab='Env1', ylab='P', type='n', las=1)

	for(i in 1:N_S){
		xvals = seq(grad[1,1],grad[2,1], length.out=100)
		yvals = sapply(xvals, function(x) niche_func(x, niches[i,'mu',1], niches[i,'sigma',1]))		
		lines(xvals, yvals, lwd=2, col=cols[i])
	}

	if(!is.null(add_env)) points(add_env[,1], rep(0, nrow(add_env)), pch=3, col=2)

	# Make plot for env var 2
	plot(c(0,0), xlim=grad[,2], ylim=c(0,1), xlab='Env2', ylab='P', type='n', las=1)

	for(i in 1:dim(niches)[1]){
		xvals = seq(grad[1,2],grad[2,2], length.out=100)
		yvals = sapply(xvals, function(x) niche_func(x, niches[i,'mu',2], niches[i,'sigma',2]))	
		lines(xvals, yvals, lwd=2, col=cols[i])
	}

	if(!is.null(add_env)) points(add_env[,2], rep(0, nrow(add_env)), pch=3, col=2)
}

# A function that plots the mutualist interaction network
#	topo = a matrix indicating the strength of interaction between mutualists
#	orderby = 'degree': ordered from most to least connected, 'name': ordered numerically by name
plot_topo = function(topo, orderby='degree', use_col='#000000', lwd=2){
	if(orderby=='degree'){
		# Calculate degree of each partner
		deg_a = rowSums(topo)
		deg_b = colSums(topo)
		
		# Order partners by degree
		ord_a = order(deg_a, decreasing=T)
		ord_b = order(deg_b, decreasing=T)
	}

	if(orderby=='name'){
		ord_a = 1:nrow(topo)
		ord_b = 1:ncol(topo)
	}
	
	# Define node locations
	pts_a = expand.grid(-1, -1*ord_a)
	pts_b = expand.grid(1, -1*ord_b)

	# Set up plot
	plot(rbind(pts_a, pts_b), type='n', axes=F, xlab='', ylab='', xlim=c(-1.5, 1.5), ylim=c(-max(dim(topo))-.5, 1))

	# Add labels
	text(-1, 0, labels='A')
	text(pts_a, labels=1:nrow(topo), pos=2)
	text(1, 0, labels='B')
	text(pts_b, labels=1:ncol(topo), pos=4)
	
	# Plot all links between pairs of potential partners
	for(i in 1:nrow(topo)){
	for(j in 1:ncol(topo)){
		shade = format(as.hexmode(floor(topo[i,j]*255)), width=2)
		segments(pts_a[i,1], pts_a[i,2], pts_b[j,1], pts_b[j,2], col=paste(use_col, shade, sep=''), lwd=lwd)
	}}

}

# A function that caculates the expected number of colonizing species for an environmental landscape.
#	niches = an S x 2 x 2 array giving mu and sigma parameters of the gaussian distribution for 2 environmental variables and S species
#	gsad = a S length vector of species global abundances
# 	grad = a list of length 2 giving the sequence of values along each gradient to be evaluated
#	draw_plot = logical indicating whether an image plot of the result should be plotted

calc_colonists = function(niches, gsad, grad = list(seq(-2,2,.2), seq(-2,2,.2)), draw_plot=F, add_sites=NA){
	N_S = dim(niches)[1]
	
	# Make array for holding probability of each species colonizing each combination of environments
	probs = array(0, dim = c(N_S, length(grad[[1]]), length(grad[[2]])))	

	# Calculate probabilities of colonization
	for(i in 1:length(grad[[1]])){
	for(j in 1:length(grad[[2]])){
		# Set environment
		env = c(grad[[1]][i], grad[[2]][j])
		
		# Calculate probability of etsablishing based on niche
		p = sapply(1:N_S, function(k) niche_func(env, niches[k,'mu',], niches[k,'sigma',]))

		# Probability of at least one propagule establishing with propagule pressure equal to global abundance
		probs[,i,j] = 1 - dbinom(0, gsad, p) 
	}}
	
	# Calculate expected number of colonist species for each time step	
	# Expected value of sum of Bernouli random variable is sum of probabilities
	ncolonists = apply(probs, c(2,3), sum)

	if(draw_plot){
		print(
		levelplot(ncolonists, xlab='Env 1', ylab='Env 2', at=seq(0,N_S, length.out=11), 
			col.regions=colorRampPalette(c('white','black'))(10), colorkey=list(labels=list(at=seq(0,N_S, length.out=11))),
			row.values = grad[[1]], column.values = grad[[2]]
		)
		)
	}
	ncolonists
}

# A function that plots a vertical color ramp on the side of a plot
# cols    : the colors to use
# n       : number of divisions
# barends : location of whole bar c(xleft, ybottom, xright, ytop)
# labels    : vector of labels for bar, assumes 1st and last numbers correspond to 1st and last colors
# title   : title to print above bar
plotColorRamp = function(cols, n, barends, labels=NA, title=NA, mycex=1.5){
	dX = barends[3] - barends[1]
	dY = barends[4] - barends[2]
	dy = dY/n
	
	xpd.old = par('xpd')
	par(xpd=T)

	usecols = colorRampPalette(cols)(n)

	for(i in 1:n){
		rect(barends[1], barends[2]+dy*(i-1), barends[3], barends[2]+dy*i, col=usecols[i], border=NA)
	}

	if(!is.na(labels)){
		dZ = labels[length(labels)]-labels[1]
		dz = dY/dZ
		Yposition = barends[2] + dz*(labels-labels[1])

		text(barends[3]+dX*0.5, Yposition, round(labels,2), pos=4, cex=mycex)
		
		segments(barends[3], Yposition, barends[3]+dX*0.5, Yposition)	
	}
	if(!is.na(title)){
		labels.round = round(labels, 2)
		
		## Determine how many characters away to place title
		digits = max(nchar(round(labels, 2))) # Maximum number of digits in a label
		largest = labels.round[which(nchar(labels.round)==digits)] # Which labels are longest
		no.decimal = sum(largest == floor(largest))>0 # Does one largest label lack a decimal?
			if(!no.decimal) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		no.negative = sum(largest >= 0)>0 # Does one largest label lack a negative sign?
			if(!no.negative) digits = digits-0.6 # Discount the size of the largest label by 0.6 a character
		
		text(barends[3]+dX*0.5+par('cxy')[1]*mycex*(digits+.5), barends[2]+0.5*dY, labels=title, srt=-90, cex=mycex)
	}
	par(xpd=xpd.old)
}






