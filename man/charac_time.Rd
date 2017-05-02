\name{charac_time}
\alias{charac_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Characteristic time measurement
}
\description{
Calculate the characteristic time it takes for the BBM+V process to reach its stationary distribution.
}
\usage{
charac_time(Npts, model)
}
\arguments{
  \item{Npts}{
The number of points used in the discretization procedure.
}
  \item{model}{
A BBM+V model fit, as returned by \link{fit_BBMV}.
}
}
\value{
The function returns the characteristic time of the process as a numeric value.
}

\author{
F. C. Boucher
}
\examples{
# Simulate data: tree + continuous trait
library(geiger)
tree=sim.bdtree(stop='taxa',n=10) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree 
# Simulate trait evolving on the tree with a linear trend towards small values 
TRAIT= Sim_BBMV(tree,x0=0,V=seq(from=0,to=5,length.out=50),sigma=10,bounds=c(-5, 5)) 
\dontrun{# fit a model with a linear  potential: 
BBM_x=fit_BBMV(tree,TRAIT,Npts=10,method='Nelder-Mead',verbose=TRUE,V_shape='linear')
# measure time to reach stationarity
charac_time(Npts=10, BBM_x)}

# compare it with tree depth
max(branching.times(tree))
}
