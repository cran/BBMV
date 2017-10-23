\name{get.landscape.FPK}
\alias{get.landscape.FPK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot Plot macroevolutionary landscapes estimated by the FPK or BBM+V models
}
\description{
Plot the adaptive landscape estimated in a BBM+V model.
}
\usage{
get.landscape.FPK(fit, Npts = 100, main = "Macroevolutionary landscape"
  , ylab = "N.exp(-V)", xlab = "Trait", xlim = NULL, ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
An FPK model fit, as returned by \link{find.mle_FPK}.
}
  \item{Npts}{
The number of points used to discretize the trait interval for plotting.
}
  \item{main}{
Title for the plot.
}
  \item{ylim}{
The upper limit of the plotting region when multiple adaptive landscapes are plotted together.
}
  \item{xlim}{
The limits of thex-axis.
}
  \item{ylab}{
Label of the y-axis.}
  \item{xlab}{
Label of the x-axis.}
}
\value{
A plot of the adaptive landscape across the trait interval.
}
\author{
F. C. Boucher
}


\examples{
\dontrun{
# Simulate data: tree + continuous trait
library(geiger)
tree=sim.bdtree(stop='taxa',n=10) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree 
# Simulate trait evolving on a macroevolutionary landscape with two peaks of equal heights
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x)) # the bounds we use for simulating
V6=10*(x^4-0.5*(x^2)+0.*x) # this is the evolutionary potential: it has two wells
TRAIT= Sim_FPK(tree,x0=0,V=V6,sigma=10,bounds=c(-5, 5)) 
# fit the FPK model: 
ll_FPK4=lnL_FPK(tree,TRAIT,Npts=25,a=NULL,b=NULL,c=NULL) # the full model
fit4=find.mle_FPK(model=ll_FPK4)
# Plot the landscape estimated
get.landscape.FPK(fit=fit4)
}

}

