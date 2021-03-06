% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_probs.r
\name{calc_probs}
\alias{calc_probs}
\title{Calculate transition probability matrices}
\usage{
calc_probs(sites, niches_a, niches_b, topo_names, poolA, poolB, mort_rate,
  assoc_probs, omega)
}
\arguments{
\item{sites}{(required) a 2 column matrix giving environmental conditions
at each site}

\item{niches_a}{(required) array of niche optima and breadths for host 
species on each environmental axis (\code{[species, optima/breadth,
environment]}) as returned by \code{\link{make_niches_gsad}}.}

\item{niches_b}{(required) array of niche optima and breadths for symbiont 
species on each environmental axis (\code{[species, optima/breadth,
environment]}) as returned by \code{\link{make_niches_gsad}}.}

\item{topo_names}{(required) matrix representing association network 
where integers label each association, as produced by 
\code{\link{name_topo}}.}

\item{poolA}{(required) binary site x species matrix indicating which host 
species are present in each site's species pool, but which are not 
necessarily established.}

\item{poolB}{(required) binary site x species matrix indicating which symbiont 
species are present in each site's species pool, but which are not 
necessarily established.}

\item{mort_rate}{(required) probability that an association dies}

\item{assoc_probs}{(required) an site x association matrix giving the 
probability that an association forms at each site if both partners 
are present. Has dimensions \code{[sites, associations]}.}

\item{omega}{(required) number controling the strength of mutualism.
\code{1-omega} is the probability that a host will establish without 
its symbiont.}
}
\value{
array of transition matrices, one for each site. 
	Has dimension \code{[association, association, site]}.
}
\description{
Calculates the probabilities that microsites will change states
for each site in a metacommunity.
}
\details{
For each site in a metacommunity, the function calculates a transition
matrix containing the probability that a microsite will changes states
during a single timestep. Microsites can be empty (0) or contain an
established host and symbiont (an 'association', given as an integer, 
\code{i}). For each potential association (given in \code{topo_names}), 
the function determines whether both host and symbiont are present in 
the species pools (\code{poolA} and \code{poolB}). If the host is not 
present, the association cannot establish (\code{P(0 -> i) = 0}). If both 
mutualists are present then the probability of an association forming is
given by \code{assoc_probs}, whereas if the symbiont is not present then
this probability is multiplied by \code{1 - omega}, which allows for 
facultative mutualism. Then, the probability that no association establishes
(\code{P(0 -> 0)}) is calculated as the product of the probabilities of  
each potential association failing to form. The remaining probability
is divided proportionally among each of the previously calculated
association probabilities (since the rows of a transition matrix must sum
to \code{1}). The probability that an association dies (\code{P(i -> 0)})
is a constant that is the same for all associations (\code{mort_rate}) 
and there is no displacement of one association by another
(\code{P(i -> j) = 0}). However, the function could be modified to allow for
competitive displacement.
}

