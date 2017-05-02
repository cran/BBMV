\name{LogLik_bounds}
\alias{LogLik_bounds}
\alias{LogLik_bounds_est}
\alias{LogLik_bounds_est_root}
\alias{bBM_loglik_0_flex_points}
\alias{bBM_loglik_bounds}
\alias{bBM_loglik_linear_bounds}
\alias{bBM_loglik_quadra_bounds}
\alias{bBM_loglik_x_flex_points}
\alias{bBM_loglik_x2x_flex_points}
\alias{bBM_loglik_x4x2x_bounds}
\alias{bBM_loglik_x4x2x_flex_pts}

\title{
Likelihood calculations for the BBM+V model
}
\description{
Internal functions that calculate the log-likelihood of various versions of the BBM+V model, used in ML and MCMC estimation. }

\usage{
LogLik_bounds(tree_formatted, dCoeff, dMat, bounds)
LogLik_bounds_est(tree, trait, dCoeff, V, bounds)
LogLik_bounds_est_root(tree, trait, dCoeff, V, x0_pos, bounds)
bBM_loglik_0_flex_points(tree, trait, Npts)
bBM_loglik_bounds(tree_formatted, dMat, bounds)
bBM_loglik_linear_bounds(tree_formatted, Npts = 100, bounds)
bBM_loglik_quadra_bounds(tree_formatted, Npts = 100, bounds)
bBM_loglik_x_flex_points(tree, trait, Npts)
bBM_loglik_x2x_flex_points(tree, trait, Npts)
bBM_loglik_x4x2x_bounds(tree_formatted, Npts = 100, bounds)
bBM_loglik_x4x2x_flex_pts(tree, trait, Npts)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree_formatted}{
A formatted tree as returned by \link{FormatTree_bounds}.
}
  \item{dCoeff}{
The diffusion coefficient.
}
  \item{dMat}{
The discretized diffusion matrix.
}
  \item{V}{
A numeric vector giving the value of the evolutionary potential in each point of the trait grid.
}
  \item{bounds}{
A vector giving the bounds of the trait interval.
}
  \item{x0_pos}{
The value of the trait at the root of the tree.
}
  \item{tree}{
A phylogenetic tree in phylo format.
}
  \item{trait}{
A trait vector for tip taxa.
}
  \item{Npts}{
The number of points used for discretizing the trait interval.
}
}

\author{
F. C. Boucher
}