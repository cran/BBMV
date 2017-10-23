\name{Uncertainty_FPK}
\alias{Uncertainty_FPK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter uncertainty}
\description{
This function plots likelihood profiles around the MLEs of paramaters and returns 95\% confidence intervals.}
\usage{
Uncertainty_FPK(fit, tree, trait, Npts = 50, effort_uncertainty = 100, 
  scope_a = c(-10, 10), scope_b = c(-10, 10), scope_c = c(-10, 10))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
An FPK model fit, as returned by \link{find.mle_FPK}.
}
  \item{tree}{
The phylogenetic tree.}
  \item{trait}{
The named trait vector}
  \item{Npts}{
The number of points used to discretize the trait interval.
}
  \item{effort_uncertainty}{
Determines the number of values at which the likelihood should be calculated for each parameter.}
  \item{scope_a}{
Extreme values that should be investigated for parameter a.}
  \item{scope_b}{
Extreme values that should be investigated for parameter b.}
  \item{scope_c}{
Extreme values that should be investigated for parameter c.}
}

\value{
A list with 95\% confidence intervals for all parameters.
}
\author{
F. C. Boucher}

\examples{
\dontrun{
# Simulate data: tree + continuous trait
library(geiger)
tree=sim.bdtree(stop='taxa',n=10) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree 
# Simulate trait evolving on a macroevolutionary landscape with two peaks of equal heights
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x)) # the bounds we use for simulating: for technical purposes only
V6=10*(x^4-0.5*(x^2)+0.*x) # this is the evolutionary potential: it has two wells
TRAIT= Sim_FPK(tree,x0=0,V=V6,sigma=10,bounds=c(-5, 5)) 
# fit the FPK model: 
ll_FPK4=lnL_FPK(tree,TRAIT,Npts=25,a=NULL,b=NULL,c=NULL) # the full model
fit4=find.mle_FPK(model=ll_FPK4)
# Measure uncertainty on model parameters
Uncertainty_FPK(fit=fit4,tree,trait=TRAIT,Npts=25,effort_uncertainty= 100,
  scope_a=c(-1,10),scope_b=c(-5,5),scope_c=c(-2,2))
}
}

