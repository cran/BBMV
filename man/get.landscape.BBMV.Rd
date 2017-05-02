\name{get.landscape.BBMV}
\alias{get.landscape.BBMV}
\alias{get.multiple.landscapes.BBMV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot macroevolutionary landscapes
}
\description{
Plot the adaptive landscape estimated in a BBM+V model.
}
\usage{
get.landscape.BBMV(model, Npts = 50, main = "Macroevolutionary landscape", 
  ylab = "-V", xlab = "Trait", xlim = NULL, ylim = NULL)
get.multiple.landscapes.BBMV(models, Npts = 50, xlim = NULL, 
  ylim = NULL, main ="Macroevolutionary landscapes", ylab = "-V", xlab ="Trait")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
A BBM+V model fit, as returned by \link{fit_BBMV}.
}
  \item{models}{
A list of BBM+V model fits, as returned by \link{fit_BBMV}.
}
  \item{landscape}{
If TRUE (the default), plots the adaptive landscape. If FALSE, plots the evolutionary potential instead.
}
  \item{Npts}{
The number of points used to discretize the trait interval for plotting.
}
  \item{main}{
A title for the plot.
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
A plot of the adaptive landscape (or alternatively the evolutionary potential) across the trait interval.
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
BBM_x=fit_BBMV(tree,TRAIT,Npts=10,method='Nelder-Mead',
  verbose=TRUE,V_shape='linear')
# plot the macroevolutionary landscape
get.landscape.BBMV(model=BBM_x,Npts=100)
}
}
