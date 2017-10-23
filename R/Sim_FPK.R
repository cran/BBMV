Sim_FPK <-
function(tree,x0=0,V=rep(0,100),sigma,bounds){
	dCoeff=log((sigma)^2/2) # the coefficient of diffusion of the model
	dMat= DiffMat_forward(V) # the transition matrix describing the probablity of evolving between two sites in the trait grid in an infinitesimal time step.
	Npts=length(V)
	ntips=length(tree$tip.label)
	trait=rep(NA,2*ntips-1) ; names(trait)=1:(2*ntips-1)
	trait[ntips+1]=x0  # root
	pMat=prep_mat_exp(dCoeff=dCoeff,dMat,bounds) # edited
for (i in 1:length(tree$edge.length)){
	proptemp= ConvProp_bounds(X= VectorPos_bounds(trait[tree$edge[i,1]],V,bounds),t=tree$edge.length[i],prep_mat = pMat) # propagate the trait forward in time: EDITED
	trait[tree$edge[i,2]]=sample(x=seq(from=bounds[1],to=bounds[2],length.out=Npts),size=1,prob= proptemp/sum(proptemp))	# sample from this probability distribution	to get a descendent node value
}
TRAIT=trait[1:ntips] ; names(TRAIT)=tree$tip.label
	return(TRAIT)
}
