\name{LogLik_bounds}
\alias{LogLik_bounds}
\alias{LogLik_bounds_est_root}
\title{
Likelihood of the FPK model
}
\description{
Internal functions use to calculate the likelihood of the FPK model, used in ML and MCMC estimation.
}
\usage{
LogLik_bounds(tree_formatted, dCoeff, dMat, bounds)
LogLik_bounds_est_root(tree, trait, dCoeff, V, x0_pos, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree_formatted}{
A formatted tree as returned by \link{FormatTree_bounds}.
}
  \item{dCoeff}{
The diffusion coefficient.
}
  \item{dMat}{
The discretized diffusion matrix.
}
  \item{V}{
A numeric vector giving the value of the evolutionary potential in each point of the trait grid.
}
  \item{bounds}{
A vector giving the bounds of the trait interval.
}
  \item{x0_pos}{
The value of the trait at the root of the tree.
}
  \item{tree}{
A phylogenetic tree in phylo format.
}
  \item{trait}{
A trait vector for tip taxa.
}
}

\author{
F. C. Boucher
}