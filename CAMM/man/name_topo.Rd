% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/name_topo.r
\name{name_topo}
\alias{name_topo}
\title{Identifies links in an association network}
\usage{
name_topo(topo)
}
\arguments{
\item{topo}{(required) association network supplied as a binary matrix}
}
\value{
matrix with the same structure as \code{topo} in which all 1s
have been replaced by integers
}
\description{
Labels each link in an association network with a unique integer.
}

