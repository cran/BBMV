\name{log_prior_nclades_plus_3_pars}
\alias{log_prior_nclades_plus_3_pars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prior function.
}
\description{
Internal function that calculates the log prior, used in MCMC estimation of the BBM+V model.
}
\usage{
log_prior_nclades_plus_3_pars(type = NULL, shape = NULL, pars, n_clades)
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
  \item{n_clades}{
The number of clades included in the multiclade analysis.}
}
\author{
F. C. Boucher
}
