\name{fit_FPK_multiple_clades_different_V_different_sig2}
\alias{fit_FPK_multiple_clades_different_V_different_sig2}
\alias{fit_BBMV_multiple_clades_different_V_different_sig2}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit independent models in several clades.}
\description{
This function is a wrapper for functions that fit the FPK or BBMV model to a single clade that does it repetitively over several clades.}
\usage{
fit_FPK_multiple_clades_different_V_different_sig2(trees, traits, 
  a = NULL, b = NULL, c = NULL, Npts = 50, method = "Nelder-Mead", init.optim = NULL)
fit_BBMV_multiple_clades_different_V_different_sig2(trees, traits, bounds, 
  a = NULL, b = NULL, c = NULL, Npts = 50, method = "Nelder-Mead", init.optim = NULL)
}
\arguments{
  \item{trees}{
A list of phylogenetic trees in 'phylo' format, one per clade.}
  \item{traits}{
A list of trait vectors for species in each clade. Should be in the same order as trees.}
 \item{bounds}{
The two bounds that constrain trait values when fitting the BBMV model. Specified by a numeric vector containing the minimum and maximum bound of the trait interval as the first and second element, respectively.
}
  \item{a}{
The value of the x^4 term in the evolutionary potential. If set to NULL (the default), this parameter will be estimated. If a numeric value is provided, this parameter will be fixed to the value specified.
}
  \item{b}{
The value of the quadratic term in the evolutionary potential. If set to NULL (the default), this parameter will be estimated. If a numeric value is provided, this parameter will be fixed to the value specified.}
  \item{c}{
The value of the linear term in the evolutionary potential. If set to NULL (the default), this parameter will be estimated. If a numeric value is provided, this parameter will be fixed to the value specified.}
  \item{Npts}{
The number of points used in the discretization procedure.
}
  \item{method}{
The optimization routine to be used: can be either "Nelder-Mead" (the default) or "L-BFGS-B". See the documentation of the \link[stats]{optim} function for more details. From our experience, "Nelder-Mead" seems to produce better results.
}
  \item{init.optim}{
A vector of initial values for model parameters to start the optimization algorithm. If left NULL (as is by default), the function chooses a reasonable starting point, but you might want to play around with it.}
}
\author{
F. C. Boucher
}

\examples{
\dontrun{
# We first create a potential that we will use to simulate trait evolution
# It has two peaks of very unequal heights
x=seq(from=-1.5,to=1.5,length.out=100)
bounds=c(min(x),max(x)) # the bounds we use for simulating
a=8 ; b=-4 ; c=1
V6=a*x^4+b*(x^2)+c*x
step_size=(max(bounds)-min(bounds))/(100-1)
V6_norm=exp(-V6)/sum(exp(-V6)*step_size) # the step size on the grid
par(mfrow=c(1,1))
plot(V6_norm,type='l')

# Now we simulate a tree and a continuous trait for 3 independent clades.
tree=sim.bdtree(stop='taxa',n=25) 
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

# Fit the FPK model using ML
fitmFPK4=fit_FPK_multiple_clades_different_V_different_sig2(trees=TREES,
  traits=TRAITS,a=NULL,b=NULL,c=NULL,Npts=50)

  }
}
