\name{Uncertainty_BBMV}
\alias{Uncertainty_BBMV}
\title{
Measure uncertainty in paramater estimation
}
\description{
This function measures the uncertainty associated with ML paramater estimates, returning plots as well as confidence intervals
}
\usage{
Uncertainty_BBMV(model, tree, trait, Npts = 50, effort_uncertainty = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A BBM+V model fit, as returned by \link{fit_BBMV}.
}
  \item{tree}{
A phylogenetic tree, in phylo format.
}
  \item{trait}{
A named vector of trait values for the tips of the tree. It should match tip labels in the phylogeny.
}
  \item{Npts}{
The number of points used to discretize the trait interval.
}
  \item{effort_uncertainty}{
The number of values of each parameter that should be tried to assess uncertainty.
}
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
\dontrun{## fit a model with a linear  potential: 
BBM_x=fit_BBMV(tree,TRAIT,Npts=10,method='Nelder-Mead',verbose=TRUE,V_shape='linear')
# Measure uncertainty around ML estimates of parameters
# Uncertainty_BBMV(BBM_x,tree,trait= TRAIT,Npts=10,effort_uncertainty= 20)
}
}