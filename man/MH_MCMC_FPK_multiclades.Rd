\name{MH_MCMC_FPK_multiclades}
\alias{MH_MCMC_FPK_multiclades}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MCMC estimation on multiple clades}
\description{
This function estimates parameter of the FPK model on multiple clades at once, making the assumption that they share the same macroevolutionary landscape but have different rates of evolution.}
\usage{
MH_MCMC_FPK_multiclades(trees, traits, bounds, Nsteps = 5e+05, record_every = 100, 
    plot_every = 500, Npts = 50, pars_init = NULL, prob_update = NULL, verbose = TRUE, 
    plot = TRUE, save_to = "MCMC_FPK_test.Rdata", save_every = 10000, type_priors = NULL, 
    shape_priors = NULL, proposal_type = "Normal", proposal_sensitivity = NULL, 
    prior.only = F, burnin.plot = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{
A list of phylogenetic trees in 'phylo' format, one per clade.}
  \item{traits}{
A list of trait vectors for species in each clade. Should be in the same order as trees.}
 \item{bounds}{
The two bounds that constrain trait values when fitting the BBMV model. Specified by a numeric vector containing the minimum and maximum bound of the trait interval as the first and second element, respectively.
}
  \item{Nsteps}{
The number of steps in the MCMC chain.}
  \item{record_every}{
How often to record a generation in the chain.}
  \item{plot_every}{
How often to plot the trace of the chain.}
  \item{Npts}{
The number of points used in the discretization procedure.
}
  \item{pars_init}{
The initial parameter values.}
  \item{prob_update}{
A numeric vector with the relative probability of update of each parameter of the model.}
  \item{verbose}{
If TRUE, prints generations to the screen.}
  \item{plot}{
If TRUE, plots the trace of parameter values along iterations during the MCMC run.}
  \item{save_to}{
The directory in which the chain should be saved.}
  \item{save_every}{
How often to save the chain.}
  \item{type_priors}{
The type of priors used, can be either normal (preferred) or uniform for log(sig2/2), a, b and c, ; and can only be discrete uniform for x0.}
  \item{shape_priors}{
A list that gives the shape for each prior. (mean,sd) for normal priors and (min,max) for continuous uniform priors. The shape is not specified for the root prior, since it is fixed to be discrete uniform on the grid.
}
  \item{proposal_type}{
The type of proposal function, only uniform is available so far.
}
  \item{proposal_sensitivity}{
The width of the uniform proposal. The entire value for x0 gives how many steps at a time can be travelled on the trait grid (better to keep it to 1)
}
  \item{prior.only}{
If TRUE, only the prior is explored but the likelihood is ignored. Default to false for estimation of the posterior.}
  \item{burnin.plot}{
The frequency of samples that should be discarded as burnin in trace plots.}
}
\details{
When specifying intial parameters yourself, be careful since x0 is actually the index of the point on the grid (between 1 and Npts_int), not the actual root value. Also the first n parameter, n being the number of clades studied, are diffusion coefficients (log(sig^2/2)), not evolutionary rates (sig^2). Finally, be careful that the bounds you propose must contain all trait values in you dataset.

Priors can be either 'Normal' (preferred) or 'Uniform' for log(sig^2/2), a, b and c. The only option for x0 is a discrete uniform prior, specified by 'Uniform'.

Each element of the shape_priors list should be a vector giving c(mean,sd) for normal priors and c(min,max) for continuous uniform priors. The shape is not specified for the root prior (it is set as 'NA' by default), since it is fixed to be discrete uniform on the grid.

Elements of the proposal_sensitivity vector can be any positive number for continuously varying parameters: c(log(sig^2/2),a,b,c). Default values should often be a good start. Only integer numbers are possible for x0 and give how many steps at a time can be travelled on the trait grid when updating these parameters. It is recommended to keep it to 1, as it is by default.}
\value{
A matrix of numeric values giving values of all parameters, the likelihood, prior and posterior at each generation sampled in the MCMC chain (one row per sample taken). The matrix has the following columns:
  \item{step }{The number of the generation sampled.}
  \item{sigsq_clade_i }{The evolutionary rate, one column per clade.}
  \item{a }{The coefficient of the x^4 term of the evolutionary potential.}
  \item{b }{The coefficient of the x^2 term of the evolutionary potential.}
  \item{c }{The coefficient of the x term of the evolutionary potential.}
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
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x)) # the bounds we use for simulating
a=8 ; b=-4 ; c=1
V6=a*x^4+b*(x^2)+c*x
step_size=(max(bounds)-min(bounds))/(100-1)
V6_norm=exp(-V6)/sum(exp(-V6)*step_size)
par(mfrow=c(1,1))
plot(V6_norm,type='l')

# Now we simulate a tree and a continuous trait for 3 independent clades. 
tree=sim.bdtree(stop='taxa',n=25) # tree with few tips for quick tests
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=1,bounds=bounds)
tree1=tree ; TRAIT1=TRAIT

tree=sim.bdtree(stop='taxa',n=25)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=0.5,bounds=bounds) 
tree2=tree ; TRAIT2=TRAIT

tree=sim.bdtree(stop='taxa',n=25)
tree$edge.length=100*tree$edge.length/max(branching.times(tree))
TRAIT= Sim_FPK(tree,x0=0.5,V=V6,sigma=0.1,bounds=bounds) 
tree3=tree ; TRAIT3=TRAIT
rm(tree) ; rm(TRAIT)

TREES=list(tree1,tree2,tree3)
TRAITS=list(TRAIT1,TRAIT2,TRAIT3)

# Fit the FPK model using ML: 
# In all clades the macroevolutionary landscape is the same 
#but they have different evolutionary rates
testbFPK4=lnl_FPK_multiclades_same_V_different_sig2(trees=TREES,
  traits=TRAITS,a=NULL,b=NULL,c=NULL,Npts=50)
fitbFPK4=find.mle_FPK_multiple_clades_same_V_different_sig2(model=testbFPK4,
  method='Nelder-Mead',init.optim=NULL)

# And now MCMC run
mcmc1=MH_MCMC_FPK_multiclades(trees=TREES,traits=TRAITS,
  bounds=fitmFPK4_SE$fits$fit_clade_1$par_fixed$bounds,Nsteps=10000,record_every=100,
  plot_every=200,Npts=25,pars_init=NULL,prob_update=NULL,verbose=TRUE,plot=TRUE,
  save_to='MCMC_FPK_test.Rdata',save_every=1000,type_priors=NULL,shape_priors=NULL,
  proposal_type='Normal',proposal_sensitivity=NULL,prior.only=F,burnin.plot=0.1)

get.landscape.FPK.MCMC(chain=mcmc1,bounds,Npts=100,burnin=0.1,probs.CI=c(0.25,0.75),
  COLOR_MEDIAN='red',COLOR_FILL='red',transparency=0.3,main='Macroevolutionary landscapes MCMC',
  ylab='N.exp(-V)',xlab='Trait',xlim=NULL,ylim=c(0,2))
  }
}
