\name{FPK_sim_traitgram}
\alias{FPK_sim_traitgram}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulations with traitgram
}
\description{
Simulates a trait evolving according to the FPK model and also plots evolution along branches of this tree using a traitgram.
}
\usage{
FPK_sim_traitgram(tree, x0, a, b, c, bounds, sigsq, time_step, res.x = 200
  , ylim.plot = NULL, return.trait = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
The phylogenetic tree on which to simulate the trait evolving.}
  \item{x0}{
The initial value of the trait (at the root of the tree).
}
  \item{a}{
The coefficient for the x^4 term of the potential.
}
  \item{b}{
The coefficient for the quadratic term of the potential.
}
  \item{c}{
The coefficient for the linear term of the potential.
}
  \item{bounds}{
The bounds on the trait interval.
}
  \item{sigsq}{
The evolutionary rate.
}
  \item{time_step}{
The time step for incremental simulations (should be smaller that the shortest branch in the tree).}
  \item{res.x}{
The number of points to use for discretizing the trait interval.
}
  \item{ylim.plot}{
The y limits of the plot. If left to NULL, the bounds of the trait interval will be used, but you might want to zoom in a bit more if the bounds are not reached.
}
  \item{return.trait}{
If set to TRUE, the function returns a named vector of trait values at the tips of the tree. }
}
\value{
The function is mainly designed to produce a plot, but eventually returns a named vector of trait values at the tips of the tree.
}
\details{
The function is slower than \link{Sim_FPK} since it simulates step by step along branches. It should be used to visualize trait evolution, but \link{Sim_FPK} should be preferred for quick simulations. 
}
\author{
F. C. Boucher}

\seealso{
\code{\link{Sim_FPK}}
}