Optim_bBM_bounds_fixed_potential <-
function(tree,trait,V,bounds=NULL,init.optim=NULL){
  if (is.null(init.optim)){
    init.optim=c(log(var(trait)/(2*max(branching.times(tree)))))
  }
  if (is.null(bounds)){
    bounds=c(min(trait),max(trait))
  }
  Npts=length(V)
  tree_formatted= FormatTree_bounds(tree,trait,V,bounds)
  dMat=DiffMat_backwards(V)
  fun= bBM_loglik_bounds(tree_formatted,dMat,bounds)
  opt=optim(par=init.optim,fn=fun,method='Brent',lower=-30,upper=10,hessian=FALSE)
  # dCoeff is log(sigma)
  # now retrieve the ML value of x0, using the ML of dCoeff
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff=opt$par,dMat,bounds)
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat) # edited
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
  }   
  x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
  ACE=lapply(tree_formatted2$Pos,function(x){cbind(seq(from=bounds[1],to=bounds[2],length.out=length(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])),x)})
  res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par),root_value=x0),lnL=-opt$value,k=4,aic=2*(4+opt$value),aicc=2*(4+opt$value)+40/(length(trait)-5),method='Brent',convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]],ACE=ACE)
  return(res)	
}
