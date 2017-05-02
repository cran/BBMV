\name{FormatTree_bounds}
\alias{FormatTree_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tree formatting.
}
\description{
Internal function used for likelihood calculation and simulation.
}
\usage{
FormatTree_bounds(tree, trait, V, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree in phylo format.
}
  \item{trait}{
A vector of traits at the tips of the tree
}
  \item{V}{
A vector giving the evolutionary potential.
}
  \item{bounds}{
A vector of bounds of the trait interval.
}
}
\author{
F. C. Boucher
}