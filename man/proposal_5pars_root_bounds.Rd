\name{proposal_5pars_root_bounds}
\alias{proposal_5pars_root_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter update for the MCMC function
}
\description{
Internal function that proposes parameter updates used in MCMC estimation of the BBMV model.
}
\usage{
proposal_5pars_root_bounds(type = "Uniform", sensitivity, pars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{type}{
The type of proposal function, only 'Uniform' is available (the default).
}
  \item{sensitivity}{
A numeric vector specifying the width of the uniform proposal for each parameter.
}
  \item{pars}{
The current parameters in the MCMC chain.
}
}
\author{
F. C. Boucher
}
