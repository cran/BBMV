\name{Optim_bBM_0_flex_pts_multiple_starts}
\alias{Optim_bBM_0_flex_pts_multiple_starts}
\alias{Optim_bBM_x4x2x}
\alias{Optim_bBM_x4x2x_flex_pts_start}
\alias{Optim_bBM_x4x2x_flex_pts_multiple_starts}
\alias{Optim_bBM_x2x_flex_pts_start}
\alias{Optim_bBM_x2x_flex_pts_multiple_starts}
\alias{Optim_bBM_x_flex_pts_start}
\alias{Optim_bBM_x_flex_pts_multiple_starts}
\alias{Optim_bBM_quadratic}
\alias{Optim_bBM_linear}
\alias{Optim_bBM_bounds_fixed_potential}
\alias{Optim_bBM_0_flex_pts_start}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum-likelihood estimation
}
\description{
Set of internal functions that are used for ML optimization of various forms of the BBM+v model.
}
\usage{
Optim_bBM_0_flex_pts_multiple_starts(tree, trait, Npts = 50, 
  method = "Nelder-Mead",verbose = T) 
Optim_bBM_x4x2x(tree, trait, Npts = 100, bounds = NULL, 
  method = "L-BFGS-B", init.optim = NULL)
Optim_bBM_x4x2x_flex_pts_start(tree, trait, Npts = 50, 
  method = "Nelder-Mead", 
	start.point=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,
	min(trait)-0.5*(max(trait) - min(trait)), 
	max(trait)+0.5*(max(trait) - min(trait))))
Optim_bBM_x4x2x_flex_pts_multiple_starts(tree, trait, Npts = 50, 
	method = "Nelder-Mead", verbose = T)
Optim_bBM_x2x_flex_pts_start(tree, trait, Npts = 50, method = "Nelder-Mead", 
	start.point = c(log(var(trait)/(2 * max(branching.times(tree)))), 0, 0, 
	min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))))
Optim_bBM_x2x_flex_pts_multiple_starts(tree, trait, Npts = 50, 
	method = "Nelder-Mead", verbose = T)
Optim_bBM_x_flex_pts_start(tree, trait, Npts = 50, method = "Nelder-Mead", 
	start.point = c(log(var(trait)/(2 * max(branching.times(tree)))), 0, 
	min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))))
Optim_bBM_x_flex_pts_multiple_starts(tree, trait, Npts = 50, 
	method = "Nelder-Mead", verbose = T)
Optim_bBM_quadratic(tree, trait, Npts = 100, bounds = NULL, 
  method = "L-BFGS-B", init.optim = NULL)
Optim_bBM_linear(tree, trait, Npts = 100, bounds = NULL, 
  method = "L-BFGS-B", init.optim = NULL)
Optim_bBM_bounds_fixed_potential(tree, trait, V, bounds = NULL, 
  init.optim = NULL)
Optim_bBM_0_flex_pts_start(tree, trait, Npts = 50, 
  method = "Nelder-Mead", 
	start.point = c(log(var(trait)/(2 * max(branching.times(tree)))), 
	min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait))))
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
  \item{bounds}{
A vector giving the bounds of the trait interval.
}
  \item{init.optim}{
A parameter combination from which to start the optimization procedure, specified as a numeric vector. Default to NULL should work in most of the cases.
}
  \item{start.point}{
A parameter combination from which to start the optimization procedure, specified as a numeric vector.}
  \item{V}{
A numeric vector giving the potential.}
}
\author{
F. C. Boucher
}