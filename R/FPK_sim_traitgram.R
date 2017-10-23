FPK_sim_traitgram <-
function(tree,x0,a,b,c,bounds,sigsq,time_step,res.x=200,ylim.plot=NULL,return.trait=FALSE){
  # reorder the tree
  tree=reorder.phylo(tree,order="cladewise") # this ordering is nicer for colors in the rainbow palette
  # initiate the plot window
  if (is.null(ylim.plot)){ylim.plot=bounds}
  COL=rainbow(dim(tree$edge)[1])
  plot(x=-max(branching.times(tree)),y=x0,main=NULL,xlim=c(-max(branching.times(tree)),0),ylim=ylim.plot,xlab='Time',ylab='Trait',pch=19,col=COL[1])
  # transform initial point and create a vector for collecting initial traits on branches to be inherited 
  init=rep(NA,(tree$Nnode+length(tree$tip.label))) ; names(init)=seq(from=1,to=(length(tree$tip.label)+tree$Nnode),by=1)
  init[length(tree$tip.label)+1]=trans_to_fixed(x0,bounds) # root trait transformed
  # get branching times of all nodes
  bt=c(rep(0,length(tree$tip.label)),-branching.times(tree)) ; names(bt)=seq(from=1,to=(length(tree$tip.label)+tree$Nnode),by=1)
  # get V and dMat
  x=seq(from=-1.5,to=1.5,length.out=res.x)
  V=a*x^4+b*x^2+c*x # potential
  dMat= DiffMat_forward(V)
  dCoeff=log(sigsq/2)
  pMat=prep_mat_exp(dCoeff,dMat,c(-1.5,1.5))
  # Now simulate along each edge of the tree
  for (i in 1:dim(tree$edge)[1]){
    n.slices=round(tree$edge.length[i]/time_step)
    temp_step=tree$edge.length[i]/n.slices # time step size is slightly modified so that we have an entire number of steps on each branch
    temp_x=rep(NA,(n.slices+1))
    temp_x[1]=init[which(names(init)==tree$edge[i,1])] # the initial trait value on this branch
    for (s in 2:length(temp_x)){ # propagate the FPK process one time step at a time
      ptemp= ConvProp_bounds(X= VectorPos_bounds(temp_x[s-1],V,c(-1.5,1.5)),t=temp_step,pMat)
      temp_x[s]=sample(x,size=1,prob=ptemp/sum(ptemp))
    }
    init[which(names(init)==tree$edge[i,2])]=temp_x[length(temp_x)]
    points(trans_from_fixed(temp_x,bounds)~seq(from=bt[which(names(bt)==tree$edge[i,1])],to=bt[which(names(bt)==tree$edge[i,2])],length.out=length(temp_x)),type='l',col=COL[i])
  }
  if(return.trait==T){
    trait=trans_from_fixed(init[1:length(tree$tip.label)],bounds) ; names(trait)=tree$tip.label
    return(trait)
  }
}
