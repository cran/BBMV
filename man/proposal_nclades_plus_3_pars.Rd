\name{proposal_nclades_plus_3_pars}
\alias{proposal_nclades_plus_3_pars}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter update for the multiclade MCMC function
}
\description{
Internal function that proposes parameter updates used in MCMC estimation of the BBMV model.
}
\usage{
proposal_nclades_plus_3_pars(type = "Uniform", sensitivity, pars, n_clades)
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
  \item{n_clades}{
The number of clades under study.
}
}
\author{
F. C. Boucher
}