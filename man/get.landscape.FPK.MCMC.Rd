\name{get.landscape.FPK.MCMC}
\alias{get.landscape.FPK.MCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot posterior distribution of macroevolutionary landscapes.}
\description{
The function plots the median value of the macroevolutionary landscape across the posterior in a solid line and draws a polygon that streches between two quantiles of the posterior.}
\usage{
get.landscape.FPK.MCMC(chain, bounds, Npts = 100, burnin = 0.1, 
  probs.CI = c(0.05, 0.95), COLOR_MEDIAN = "red", COLOR_FILL = "red", 
  transparency = 0.3, main = "Macroevolutionary landscapes MCMC", 
  ylab = "N.exp(-V)", xlab = "Trait", xlim = NULL, ylim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{
An data.frame object representing the output of an MCMC chain, as obtained by \link{MH_MCMC_FPK}.}
  \item{bounds}{
The bounds on the trait interval}
  \item{Npts}{
The number of points used in the discretization procedure.
}
  \item{burnin}{
The percentage of generations discarded as burnin.}
  \item{probs.CI}{
A vector of the two quantiles of the posterior distribution between which samples should be considered.}
  \item{COLOR_MEDIAN}{
The color used to plot the median macroevolutionary landscape across the posterior.}
  \item{COLOR_FILL}{
The color used to plot the polygon that stretches between the two quantiles of the posterior.}
  \item{transparency}{
The transparency used for plotting the polygon}
  \item{main}{
Title of the graph.}
  \item{ylab}{
y label of the graph.}
  \item{xlab}{
X label of the graph.}
  \item{xlim}{
%%     ~~Describe \code{xlim} here~~
}
  \item{ylim}{
%%     ~~Describe \code{ylim} here~~
}
}
\author{
F.C. Boucher}

\seealso{
\link{MH_MCMC_FPK}}
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
get.landscape.FPK.MCMC(chain=MCMC,bounds=c(5,5),Npts=100,burnin=0.1,
  probs.CI=c(0.025,0.975),COLOR_MEDIAN='red',COLOR_FILL='red',transparency=0.3,
  main='Macroevolutionary landscapes MCMC',ylab='N.exp(-V)',xlab='Trait',
  xlim=NULL,ylim=NULL)  
}
}

