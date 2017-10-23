\name{Sim_FPK}
\alias{Sim_FPK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulation of the BBM+V process.
}
\description{
The function simulates a continuous trait evolving according to the FPK process along the branches of a phylogenetic tree.
}
\usage{
Sim_FPK(tree, x0 = 0, V = rep(0, 100), sigma, bounds)
}
\arguments{
  \item{tree}{
A phylogenetic tree in phylo format.
}
  \item{x0}{
The value of the trait at the root of the tree.
}
  \item{V}{
A vector giving the values of the evolutionary potential at each point of the discretized trait grid. Default is a flat potential, i.e. bounded Brownian Motion.
}
  \item{sigma}{
The square root of the diffusion rate.
}
  \item{bounds}{
A vector giving the values of the bounds of the trait interval.
}
}

\value{
A numeric vector with values of the trait at the tips of the tree. Names correspond to tip labels in the tree.
}

\author{
F. C. Boucher
}

\examples{
# Simulate data: tree + continuous trait
library(geiger)
tree=sim.bdtree(stop='taxa',n=20) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree
# Simulate a trait evolving on the tree with a linear trend towards small values
TRAIT= Sim_FPK(tree,x0=0,V=seq(from=0,to=5,length.out=50),sigma=10,bounds=c(-5, 5)) 
}
