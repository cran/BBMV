\name{log_prior_7pars_root_bounds}
\alias{log_prior_7pars_root_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior function.
}
\description{
Internal function that calculates the log prior, used in MCMC estimation of the BBM+V model.
}
\usage{
log_prior_7pars_root_bounds(type = c(rep("Normal", 4), rep("Uniform", 3)), 
	shape = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10), NA, 10, 10), 
	pars, Npts_int, trait)
}
%- maybe also 'usage' for other objects documented here.
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
  \item{Npts_int}{
The number of points on the grid between min(trait) and max(trait).
}
  \item{trait}{
A named vector of trait values for the tips of the tree.
}
}

\author{
F. C. Boucher
}
