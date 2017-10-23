\name{lnL_FPK}
\alias{lnL_FPK}
\alias{lnL_BBMV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creation of the likelihood function
}
\description{
Functions that builds the likelihood function of the FPK or BBMV model
}
\usage{
lnL_FPK(tree, trait, a = NULL, b = NULL, c = NULL, Npts)
lnL_BBMV(tree, trait,bounds, a = NULL, b = NULL, c = NULL, Npts)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree in 'phylo' format
}
  \item{trait}{
A named vector of trait values for the tips of the tree. It should match tip labels in the phylogeny. Alternatively, a named list with one element per tip in the tree, each element being in turn a numeric vector with multiple measurements of the trait for this tip.}
  \item{bounds}{
The two bounds that constrain trait values when fitting the BBMV model. Specified by a numeric vector containing the minimum and maximum bound of the trait interval as the first and second element, respectively.
}
  \item{a}{
The value of the x^4 term in the evolutionary potential. If set to NULL (the default), this parameter will be estimated. If a numeric value is provided, this parameter will be fixed to the value specified.
}
  \item{b}{
The value of the quadratic term in the evolutionary potential. If set to NULL (the default), this parameter will be estimated. If a numeric value is provided, this parameter will be fixed to the value specified.}
  \item{c}{
The value of the linear term in the evolutionary potential. If set to NULL (the default), this parameter will be estimated. If a numeric value is provided, this parameter will be fixed to the value specified.}
  \item{Npts}{
The number of points used in the discretization procedure.
}
}

\value{
A list of several items, including the data and model call, but most importantly the likelihood function ($fun element).}

\author{
F.C. Boucher}

\seealso{
\code{\link{find.mle_FPK}}
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
# create a likelihood function for the FPK model: 
ll_FPK4=lnL_FPK(tree,TRAIT,Npts=25,a=NULL,b=NULL,c=NULL) # the full model
}
}