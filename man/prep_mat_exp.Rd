\name{prep_mat_exp}
\alias{prep_mat_exp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Matrix exponential.
}
\description{
Internal function used for likelihood calculation and simulation.
}
\usage{
prep_mat_exp(dCoeff, dMat, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dCoeff}{
The diffusion coefficient.
}
  \item{dMat}{
The diffusion matrix.
}
  \item{bounds}{
A vector with two bounds for the trait interval.
}
}
\author{
F. C. Boucher
}