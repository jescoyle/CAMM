% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summarize_camm.r
\name{summarize_camm}
\alias{summarize_camm}
\title{Summarize simulation results}
\usage{
summarize_camm(results, what, type = NA)
}
\arguments{
\item{results}{(required) list of simulations results, as returned by
\code{link{run_sim_N}}}

\item{what}{(required) character string indicating which community 
descriptor to summarize (see details)}

\item{type}{character string indicating which community to summarize
(see details: required if \code{what} is \code{'S'}, \code{'N'}, or 
\code{'cor'})}
}
\value{
an array
}
\description{
Summarizes multiple independent simulation runs
}
\details{
This function summarizes a given community descriptor across multiple
independent simulations run on the same set of parameters. It extracts
the requested descriptor for each simulation run in \code{results}
and then calculates the mean, variance, median and 2.5th and 97.5th 
quantiles across simulation runs.

The following community descriptors can be requested 
(using the paramter \code{what}):
\describe{
\item{'S'}{mean site species richness, total species richness, and 
	beta diversity}
\item{'N'}{mean site abundance and total abundance}
\item{'Cor_ab'}{correlation between host and symbiont richness}
\item{'cor'}{community ~ environment correlations}
}
The user must also indicate which community should be summarized
useing \code{type}: established associations (\code{'species'}), 
or host (\code{'a'}) or symbiont (\code{'b'}) species involved in
established associations.
}

