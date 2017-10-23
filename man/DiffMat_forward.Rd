\name{DiffMat_forward}
\alias{DiffMat_forward}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Diffusion matrix building
}
\description{
Internal function that builds the discretized diffusion matrix of the FPK process going forward in time (for simulations)
}
\usage{
DiffMat_forward(V)
}
\arguments{
  \item{V}{
A vector giving the values of the evolutionary potential (V) at each point in the gridded trait interval.
}
}
\author{
F.C. Boucher
}
