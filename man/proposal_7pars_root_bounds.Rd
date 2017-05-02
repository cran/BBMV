\name{proposal_7pars_root_bounds}
\alias{proposal_7pars_root_bounds}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameter update
}
\description{
Internal function that proposes parameter updates used in MCMC estimation of the BBM+V model.
}
\usage{
proposal_7pars_root_bounds(type = "Uniform", sensitivity, pars, trait, Npts_int)
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
  \item{trait}{
A named vector of trait values for the tips of the tree. It should match tip labels in the phylogeny.
}
  \item{Npts_int}{
The number of points on the grid between min(trait) and max(trait).
}
}
\author{
F. C. Boucher
}
