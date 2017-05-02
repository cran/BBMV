\name{fit_BBMV}
\alias{fit_BBMV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to fit the BBM+V model
}
\description{
Fits various forms of the BBM+V model to trait data and a phylogeny, using maximum likelihood.
}
\usage{
fit_BBMV(tree, trait, Npts = 50, method = "Nelder-Mead", 
  verbose = T, V_shape, bounds = "estimate", init.optim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree in phylo format
}
  \item{trait}{
A named vector of trait values for the tips of the tree. It should match tip labels in the phylogeny.
}
  \item{Npts}{
The number of points used in the discretization procedure.
}
  \item{method}{
The optimization routine used.
}
  \item{verbose}{
If TRUE, prints progress to the screen.
}
  \item{V_shape}{
The shape of the evolutionary potential. This should be one of 'flat','linear','quadratic', or 'full'. Alternatively, a given potential can be specified using a vector of values for V, along with bounds. 
}  \item{bounds}{
Describes how to treat the bounds. If set to 'estimate' (the default), bounds will be estimated along with other parameters. Alternatively, you can fix two bounds by specifying a numeric vector containing the minimum and maximum bound of the trait interval as the first and second element, respectively.}
  \item{init.optim}{
A parameter combination from which to start the optimization procedure, specified as a numeric vector. Default to NULL should work in most of the cases.}
}
\details{
Only 'Nelder-Mead' (the default, which we found to be more reliable) and 'L-BFGS-B' are supported as optimization methods. For further details on both methods see the help of the \link[stats]{optim} function.

The different values for V_shape refer to the following shapes of the evolutionary potential: (1) 'flat' is a flat potential as in BBM, i.e. \eqn{V(x)=0} ; (2) 'linear' is a linear potential, i.e. \eqn{V(x)=c.x} ; (3) 'quadratic' is a quadratic potential, i.e. \eqn{V(x)=b.x^2+c.x} ; (4) 'full' is a potential with three polynomial terms, i.e. \eqn{V(x)=a.x^4+b.x^2+c.x}.
}

\value{
A list with the following elements:
  \item{par }{The ML estimates of model parameters: the two bounds of the trait interval ($bounds), the diffusion rate ($sigsq), the coefficients of the evolutionary potential ($a, $b, and $c, depending on the shape fitted), and the value of the trait at the root of the tree ($root_value). }
  \item{lnL }{The log-likelihood of the model.}
  \item{k }{The number of parameters of the model.}
  \item{aic }{The AIC of the model.}
  \item{aicc }{The AICc of the model.}
  \item{method }{The optimization method used to fit the model.}
  \item{convergence }{Convergence code returned by optim. 0 indicates successful convergence. For other values see the help of the \link[stats]{optim} function.}
  \item{message }{Convergence message returned by optim. See the help of the \link[stats]{optim} function.}
  \item{root_density }{A vector giving the density of probability of the trait at the root of the tree. Each element corresponds to one point in the discretized trait grid, from left (lower bound) to right (upper bound).}
  \item{ACE}{A list of tables, one for each node of the phylogeny. Each table has the position of the points in the trait grid as its first column, and their associated probabilities in the second column.}
}
\author{
F. C. Boucher
}
\examples{
# Simulate data: tree + continuous trait
library(geiger)
tree=sim.bdtree(stop='taxa',n=10) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree 
# Simulate trait evolving on the tree with a linear trend towards small values 
TRAIT= Sim_BBMV(tree,x0=0,V=seq(from=0,to=5,length.out=50),sigma=10,bounds=c(-5, 5)) 
\dontrun{# fit a model with a linear  potential: 
BBM_x=fit_BBMV(tree,TRAIT,Npts=10,method='Nelder-Mead',verbose=TRUE,V_shape='linear')}
}
