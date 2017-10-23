find.mle_FPK <-
function(model,method='Nelder-Mead',init.optim=NULL,safe=F){
  # safe=T for safer optimization starting from 3 different starting points
  if (safe==F){ # only one optimization
    if(is.null(init.optim)==T){init.optim=c(-10,rep(0,model$ncoeff))}
    else{}
    opt=optim(par=init.optim,fn=model$fun,method=method,control=list(maxit=50000))
  }
  else {
    init=c(-10,-1,0) #diffusion Coefficient
    starts=list()
    for (i in 1:3){
      starts[[i]]=optim(par=c(init[i],rep(0,model$ncoeff)),fn=model$fun,method=method,control=list(maxit=50000))
      cat("Finished preliminary fit, log-lik= ",-starts[[i]]$value,sep="\n")
    }
    lnls=c(-starts[[1]]$value)
    for (i in 2:length(starts)){lnls=c(lnls,-starts[[i]]$value)}
    opt=starts[[which(lnls==max(lnls))]] # the starting point which gave the highest lnl
  }
  par_fixed=model$par_fixed
  par=list()
  par$sigsq=2*exp(opt$par[1])
  if (model$ncoeff==3){par$a=opt$par[2] ; par$b=opt$par[3] ; par$c=opt$par[4]}
  else if (model$ncoeff==2){par$b=opt$par[2] ; par$c=opt$par[3]}
  else if (model$ncoeff==1){par$c=opt$par[2]}
  else{}
  # now retrieve x0
  # retrieve coefficients for the potential
  if ('a'%in%names(par)){a=par$a}
  else {a=par_fixed$a}
  if ('b'%in%names(par)){b=par$b}
  else {b=par_fixed$b}
  if ('c'%in%names(par)){c=par$c}
  else {c=par_fixed$c}
  bounds=par_fixed$bounds
  Npts=model$Npts ; tree=model$tree ; trait=model$trait
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  V=a*SEQ^4+b*SEQ^2+c*SEQ
  dMat=DiffMat_backwards(V)
  ll_root=as.data.frame(matrix(NA, Npts,2))
  colnames(ll_root)=c('root','density')
  ll_root$root =seq(bounds[1],to=bounds[2], length.out= Npts)
  tree_formatted= FormatTree_bounds(tree,trait,V=rep(0, Npts),bounds=bounds)
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff=log(par$sigsq/2),dMat,bounds=bounds) # edited
  logFactor=0
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
    norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
    logFactor=logFactor+log(norm)
  }
  ll_root$density =tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]
  return(res=list(lnL=-opt$value,aic=2*(length(init.optim)+opt$value),k=length(init.optim),par=par,par_fixed=par_fixed,root=ll_root,convergence=opt$convergence,message=opt$message,tree=tree,trait=trait,Npts=Npts)) 
}
