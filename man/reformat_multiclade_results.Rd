\name{reformat_multiclade_results}
\alias{reformat_multiclade_results}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Format the output of a multiclade fit
}
\description{
This functions takes the output of a multiclade fit and formats it as a list of model fits for each clade. Functions used to analyze single clade fits can then be used.
}
\usage{
reformat_multiclade_results(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
A multiclade model fit, as returned by \link{find.mle_FPK_multiple_clades_same_V_different_sig2}, \link{find.mle_FPK_multiple_clades_same_V_same_sig2}, \link{fit_BBMV_multiple_clades_different_V_different_sig2} or \link{fit_FPK_multiple_clades_different_V_different_sig2}.
}
}
\value{
A list containing model fits for each clade, in the same format as the object returned by \link{find.mle_FPK}.
}
\author{
F.C. Boucher}
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

# Fit the FPK model
testbFPK4=lnl_FPK_multiclades_same_V_different_sig2(trees=TREES,
  traits=TRAITS,a=NULL,b=NULL,c=NULL,Npts=50)
fitbFPK4=find.mle_FPK_multiple_clades_same_V_different_sig2(model=testbFPK4,
  method='Nelder-Mead',init.optim=NULL)
fits=reformat_multiclade_results(fitbFPK4)
}
}

