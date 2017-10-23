\name{MH_MCMC_FPK}
\alias{MH_MCMC_FPK}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC estimation
}
\description{
The function estimates the parameters of the BBM+V model using an MCMC chain with the Metropolis Hastings algorithm and a Gibbs sampler.
}
\usage{
MH_MCMC_FPK(tree, trait, bounds, Nsteps = 5e+05, record_every = 100, plot_every = 500, 
  Npts = 50, pars_init = c(0, 0, 0, 0, 25), prob_update = c(0.2, 0.2, 0.2, 0.2, 0.2), 
  verbose = TRUE, plot = TRUE, save_to = "MCMC_FPK_test.Rdata", save_every = 10000, 
  type_priors = c(rep("Normal", 4), "Uniform"), 
  shape_priors = list(c(0, 10), c(0, 10), c(0, 10), c(0, 10), NA), 
  proposal_type = "Uniform", proposal_sensitivity = c(0.1, 0.1, 0.1, 0.1, 1), 
  prior.only = F, burnin.plot = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{
A phylogenetic tree in phylo format.
}
  \item{trait}{
A named vector of trait values for the tips of the tree. It should match tip labels in the phylogeny.
}
  \item{bounds}{
A vector with two elements giving the bounds on the trait interval.}
  \item{Nsteps}{
The number of generations in the MCMC chain.
}
  \item{record_every}{
The frequency used for sampling the MCMC chain.
}
  \item{plot_every}{
The frequency at which the chain is plotted (if plot=TRUE).
}
  \item{Npts}{
The number of points on the grid between the bounds.
}
  \item{pars_init}{
A vector giving the initial parameters for starting the algorithm, which correspond to the following: c(log(sig^2/2),a,b,c,x0).
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
  \item{burnin.plot}{
The percentage of samples considered as burnin and thus not shown on the trace plot that the function produces.}
}
\details{
When specifying intial parameters yourself, be careful since x0 is actually the index of the point on the grid (between 1 and Npts_int), not the actual root value. Also the first parameter is the diffusion coefficient (log(sig^2/2)), not the evolutionary rate (sig^2). Finally, be careful that the bounds you propose must contain all trait values in you dataset.

Priors can be either 'Normal' (preferred) or 'Uniform' for log(sig^2/2), a, b and c. The only option for x0 is a discrete uniform prior, specified by 'Uniform'.

Each element of the shape_priors list should be a vector giving c(mean,sd) for normal priors and c(min,max) for continuous uniform priors. The shape is not specified for the root prior (it is set as 'NA' by default), since it is fixed to be discrete uniform on the grid.

Elements of the proposal_sensitivity vector can be any positive number for continuously varying parameters: c(log(sig^2/2),a,b,c). Default values should often be a good start. Only integer numbers are possible for x0 and give how many steps at a time can be travelled on the trait grid when updating these parameters. It is recommended to keep it to 1, as it is by default.}
\value{
A matrix of numeric values giving values of all parameters, the likelihood, prior and posterior at each generation sampled in the MCMC chain (one row per sample taken). The matrix has the following columns:
  \item{step }{The number of the generation sampled.}
  \item{sigsq }{The evolutionary rate.}
  \item{a }{The coefficient of the x^4 term of the evolutionary potential.}
  \item{b }{The coefficient of the x^2 term of the evolutionary potential.}
  \item{c }{The coefficient of the x term of the evolutionary potential.}
  \item{root }{The value of the trait at the root of the tree.}
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
\dontrun{
# Simulate data: tree + continuous trait
library(geiger)
tree=sim.bdtree(stop='taxa',n=10) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree)) # rescale the tree 
# Simulate trait evolving on a macroevolutionary landscape with two peaks of equal heights
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x)) # the bounds we use for simulating: for technical purposes only
V6=10*(x^4-0.5*(x^2)+0.*x) # this is the evolutionary potential: it has two wells
TRAIT= Sim_FPK(tree,x0=0,V=V6,sigma=10,bounds=c(-5, 5)) 
# Run a MCMC chain to fit the FPK model
MCMC=MH_MCMC_FPK(tree,trait=TRAIT,bounds=c(5,5),Nsteps=10000,record_every=100,
  plot_every=100,Npts=20,pars_init=c(0,-4,-4,0,1),prob_update=c(0.2,0.25,0.25,0.25,0.05),
  verbose=TRUE,plot=TRUE,save_to='MCMC_FPK_test.Rdata',save_every=100,
  type_priors=c(rep('Normal',4),'Uniform'),
  shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',
  proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F)
}
}