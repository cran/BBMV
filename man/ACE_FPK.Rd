\name{ACE_FPK}
\alias{ACE_FPK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Ancestral Character Estimation
}
\description{
Function to perform Ancestral Character Estimation under the FPK (or BBM+V) model
}
\usage{
ACE_FPK(fit, specific.point = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
An FPK model fit, as returned by \link{find.mle_FPK}.
}
  \item{specific.point}{
If set to NULL (the default), then the function will produce an ACE at all internal nodes in the tree. Alternatively, specific.point can be used to ask for an ACE at any specific point in the tree (i.e., not a node): specific.point must then be a vector with three elements: c(parent_node,child_node,time_from_start_of branch).
}
}
\value{
For each internal node, the function returns a table giving the probability density of the trait. The first column gives all possible trait values on the discretized trait grid, and the second the probability density at each of these points. If only a specific.point was asked, the function only returns one such table.
}
\author{
F. C. Boucher
}