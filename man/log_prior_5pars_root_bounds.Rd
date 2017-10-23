\name{log_prior_5pars_root_bounds}
\alias{log_prior_5pars_root_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior function.
}
\description{
Internal function that calculates the log prior, used in MCMC estimation of the BBM+V model.
}
\usage{
log_prior_5pars_root_bounds(type = c(rep("Normal", 4), "Uniform") , shape = 
  list(c(0, 10), c(0, 10), c(0, 10), c(0, 10), NA), pars, Npts)
}
\arguments{
  \item{type}{
A vector giving the type of prior for each parameter.
}
  \item{shape}{
A list giving the shape of the prior for each parameter.
}
  \item{pars}{
The parameter values at which the prior should be evaluated.
}
  \item{Npts}{
The number of points on the grid.
}
}

\author{
F. C. Boucher
}
