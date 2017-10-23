\name{charac_time}
\alias{charac_time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Characteristic time measurement
}
\description{
Calculate the characteristic time it takes for the FPK process to reach its stationary distribution.
}
\usage{
charac_time(Npts, fit)
}
\arguments{
  \item{Npts}{
The number of points used in the discretization procedure.
}
  \item{fit}{
A FPK model fit, as returned by \link{find.mle_FPK}.
}
}
\value{
The function returns the characteristic time of the process as a numeric value.
}

\author{
F. C. Boucher
}

