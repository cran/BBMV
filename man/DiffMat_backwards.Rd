\name{DiffMat_backwards}
\alias{DiffMat_backwards}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diffusion matrix building
}
\description{
Internal function that builds the discretized diffusion matrix of the BBM+V process going backwards in time (for likelihood calculations)
}
\usage{
DiffMat_backwards(V)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
A vector giving the values of the evolutionary potential (V) at each point in the gridded trait interval.
}
}
\author{
F. C. Boucher
}
