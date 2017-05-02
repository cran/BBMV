\name{MH_MCMC_V_ax4bx2cx_root_bounds}
\alias{MH_MCMC_V_ax4bx2cx_root_bounds}
\title{
MCMC estimation
}
\description{
The function estimates the parameters of the BBM+V model using an MCMC chain with the Metropolis Hastings algorithm and a Gibbs sampler.
}
\usage{
MH_MCMC_V_ax4bx2cx_root_bounds(tree, trait, Nsteps = 5e+05, record_every = 100, 
	plot_every = 500, Npts_int = 50, pars_init = c(0, 0, 0, 0, 25, -10, 10), 
	prob_update = c(0.25, 0.2, 0.2, 0.2, 0.05, 0.05, 0.05), verbose = TRUE, 
	plot = TRUE, save_to = "testMCMC.Rdata", save_every = 10000, 
	type_priors = c(rep("Normal", 4), rep("Uniform", 3)), 
	shape_priors = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10), NA, 10, 10), 
	proposal_type = "Uniform", proposal_sensitivity = c(0.1, 0.1, 0.1, 0.1, 1, 1, 1), 
	prior.only = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree in phylo format.
}
  \item{trait}{
A named vector of trait values for the tips of the tree. It should match tip labels in the phylogeny.
}
  \item{Nsteps}{
The number of generations in the MCMC chain.
}
  \item{record_every}{
The frequency used for sampling the MCMC chain.
}
  \item{plot_every}{
The frequency at which the chain is plotted (if plot=TRUE).
}
  \item{Npts_int}{
The number of points on the grid between min(trait) and max(trait).
}
  \item{pars_init}{
A vector giving the initial parameters for starting the algorithm, which correspond to the following: c(log(sig^2/2),a,b,c,x0,Bmin,Bmax).
}
  \item{prob_update}{
A vector giving the relative frequencies of update of the different parameters of the model.
}
  \item{verbose}{
If TRUE, will print some generations of the chain to the screen.
}
  \item{plot}{
If TRUE, the chain is plotted from time to time.
}
  \item{save_to}{
The path to the file where the chain is saved (can be useful in case the chain crashes).
}
  \item{save_every}{
Sets how often the chain is saved.
}
  \item{type_priors}{
A character vector specifying the type of priors used. Either 'Uniform' or 'Normal'. See Details.
}
  \item{shape_priors}{
A list that gives the shape for each prior. See Details.
}
  \item{proposal_type}{
The type of proposal function, only 'Uniform' is available (the default).
}
  \item{proposal_sensitivity}{
A numeric vector specifying the width of the uniform proposal for each parameter. See Details.
}
  \item{prior.only}{
Default to FALSE for estimation of the posterior. If TRUE, the likelihood is not evaluated: this is mostly useful for internal test of the Gibbs sampler.
}
}
\details{
When specifying intial parameters yourself, be careful since x0 is actually the index of the point on the grid (between 1 and Npts_int), not the actual root value. Also the fist parameter (log(sig^2/2)) is the diffusion coefficient, not the evolutionary rate (sig^2). Finally, be careful that the intial bounds you propose actually contain all trait values in you dataset.

Priors can be either 'Normal' (preferred) or 'Uniform' for log(sig^2/2), a, b and c. The only option for for bounds and x0 is discrete uniform priors, specified by 'Uniform'.

Each element of the shape_priors list should be a vector giving c(mean,sd) for normal priors and c(min,max) for continuous uniform priors. The shape is not specified for the root prior (it is set as 'NA' by default), since it is fixed to be discrete uniform on the grid. Values for the priors on the bounds (forced to be discrete uniform) are a single numeric value giving the maximum number of points that can be added on the trait grid outside of the observed trait interval, with a default value of 10 points.

Elements of the proposal_sensitivity vector can be any positive number for continuously varying parameters: c(log(sig^2/2),a,b,c). Default values should often be a good start. Only integer numbers are possible for x0, Bmin, and Bmax and give how many steps at a time can be travelled on the trait grid when updating these parameters. It is recommended to keep it to 1, as it is by default.
}

\value{
A matrix of numeric values giving values of all parameters, the likelihood, prior and posterior at each generation sampled in the MCMC chain (one row per sample taken). The matrix has the following columns:
  \item{step }{The number of the generation sampled.}
  \item{sigsq }{The evolutionary rate.}
  \item{a }{The coefficient of the x^4 term of the evolutionary potential.}
  \item{b }{The coefficient of the x^2 term of the evolutionary potential.}
  \item{c }{The coefficient of the x term of the evolutionary potential.}
  \item{root }{The value of the trait at the root of the tree.}
  \item{bmin }{The value of the lower bound of the trait interval.}
  \item{bmax }{The value of the upper bound of the trait interval.}
  \item{lnprior }{The logarithm of the prior.}
  \item{lnlik }{The logarithm of the likelihood.}
  \item{quasi-lnpost }{The logarithm of the (unnormalized) posterior.}
  \item{Acceptance }{Whether the proposed MCMC move was accepted (1) or not (0).}
  \item{Par_updated }{Which parameter was updated in this generation.}
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
\dontrun{
# fit the full model using MCMC estimation: 
MCMC= MH_MCMC_V_ax4bx2cx_root_bounds(tree,trait=TRAIT,Nsteps=20000,
  record_every=100,plot_every=500,Npts_int=20,pars_init=c(-8,0,0,0,5,min(TRAIT),max(TRAIT)),
  prob_update=c(0.05,0.3,0.3,0.15,0.15,0.05,0.05),verbose=TRUE,plot=TRUE,
  save_to='testMCMC.Rdata',save_every=1000,type_priors=c(rep('Normal',4),rep('Uniform',3)),
  shape_priors=list(c(0,2),c(0,2),c(0,2),c(0,2),NA,30,30),proposal_type='Uniform',
  proposal_sensitivity=c(1,0.5,0.5,0.5,1,1,1),prior.only=F)}
}