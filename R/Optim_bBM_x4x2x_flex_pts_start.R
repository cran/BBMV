Optim_bBM_x4x2x_flex_pts_start <-
function(tree,trait,Npts=50,method='Nelder-Mead',start.point=c(log(var(trait)/(2*max(branching.times(tree)))),0,0,0,min(trait)-0.5*(max(trait)-min(trait)),max(trait)+0.5*(max(trait)-min(trait)))){ 	# Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
  # Nelder-Mead is the default method since it seems to be more robust with this rather high number of parameters
  if (!(method%in%c('L-BFGS-B','Nelder-Mead'))){stop('Incorrect optimization method')}
  fun= bBM_loglik_x4x2x_flex_pts(tree,trait,Npts)
  if (method=='L-BFGS-B'){
    opt=optim(par=start.point,fn=fun,method='L-BFGS-B',lower=c(-10,-10,-10,-10,min(trait)-10*(max(trait)-min(trait)),max(trait)),upper=c(10,10,10,10,min(trait),max(trait)+10*(max(trait)-min(trait))),hessian=FALSE)
  }
  if (method=='Nelder-Mead'){
    opt=optim(par=start.point,fn=fun,method='Nelder-Mead',hessian=FALSE,control=list(maxit=10000))
  }
  # dCoeff is log(sigma^2/2)
  # now retrieve the ML value of x0, using the ML of dCoeff, bounds and V
  SEQ=seq(from=0,to=1,length.out=Npts)
  V=opt$par[2]*SEQ^4+opt$par[3]*SEQ^2+opt$par[4]*SEQ
  dMat=DiffMat_backwards(V)
  bounds=c(opt$par[5],opt$par[6])
  tree_formatted=FormatTree_bounds(tree,trait,V,bounds)
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff=opt$par[1],dMat,bounds) # edited
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
  }   
  x0=bounds[1]+(bounds[2]-bounds[1])*(which(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]==max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))-1)/(Npts-1)
  ACE=lapply(tree_formatted2$Pos,function(x){cbind(seq(from=bounds[1],to=bounds[2],length.out=length(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])),x)})
  res=list(par=list(bounds=bounds,sigsq=2*exp(opt$par[1]),a=opt$par[2],b=opt$par[3],c=opt$par[4],root_value=x0),lnL=-opt$value,k=7,aic=2*(7+opt$value),aicc=2*(7+opt$value)+(112)/(length(trait)-8),method=method,convergence=opt$convergence,message=opt$message,root_density=tree_formatted2$Pos[[tree_formatted2$tab[i,1]]],ACE=ACE)
  return(res)	
  
}
