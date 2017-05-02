\name{ConvProp_bounds}
\alias{ConvProp_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convolution of the diffusion matrix with the trait density vector.
}
\description{
Internal function used for likelihood calculation and simulation.
}
\usage{
ConvProp_bounds(X, t, prep_mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
A trait density vector.
}
  \item{t}{
The time over which to do the convolution (usually the length of one branch).
}
  \item{prep_mat}{
The diagonalized diffusion matrix.
}
}
\author{
F. C. Boucher
}

