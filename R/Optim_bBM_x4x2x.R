Optim_bBM_x4x2x <-
function(tree,trait,Npts=100,bounds=NULL,method='L-BFGS-B',init.optim=NULL){
  if (is.null(init.optim)){
    init.optim=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0)
  }
  if (is.null(bounds)){
    bounds=c(min(trait),max(trait))
  }
  if (!(method%in%c('L-BFGS-B','Nelder-Mead'))){stop('Incorrect optimization method')}
  tree_formatted= FormatTree_bounds(tree,trait,rep(0,Npts),bounds) # we don't care about the shape of the potential to format tree & trait
  fun= bBM_loglik_x4x2x_bounds(tree_formatted,Npts=Npts,bounds)
  if (method=='L-BFGS-B'){
    opt=optim(par=init.optim,fn=fun,method='L-BFGS-B',lower=c(-10,-10,-10,-10),upper=c(10,10,10,10),hessian=FALSE)
  }
  if (method=='Nelder-Mead'){
    opt=optim(par=init.optim,fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
  }
  # dCoeff is log(sigma^2/2)
  # now retrieve the ML value of x0, using the ML of dCoeff
  tree_formatted2= tree_formatted
  SEQ=seq(from=0,to=1,length.out=Npts)
  V=opt$par[2]*SEQ^4+opt$par[3]*SEQ^2+opt$par[4]*SEQ
  dMat=DiffMat_backwards(V)
  pMat=prep_mat_exp(dCoeff=opt$par[1],dMat,bounds) # edited
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat) #edited
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
  }   
  x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
  ACE=lapply(tree_formatted2$Pos,function(x){cbind(seq(from=bounds[1],to=bounds[2],length.out=length(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])),x)})
  res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],b=opt$par[3],c=opt$par[4],root_value=x0),lnL=-opt$value,k=7,aic=2*(7+opt$value),aicc=2*(7+opt$value)+(112)/(length(trait)-8),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]],ACE=ACE)
  return(res)	
}
