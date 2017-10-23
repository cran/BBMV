\name{VectorPos_bounds}
\alias{VectorPos_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Discretization of a continuous trait value into a probability vector.
}
\description{
Internal function used for likelihood calculation and simulation.
}
\usage{
VectorPos_bounds(x, V, bounds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A numeric value of the trait or a vector containing multiple measurements.
}
  \item{V}{
The evolutionary potential used
}
  \item{bounds}{
A vector with the values of both bounds.
}
}
\author{
F. C. Boucher
}
