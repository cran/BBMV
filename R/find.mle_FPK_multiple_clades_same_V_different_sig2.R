find.mle_FPK_multiple_clades_same_V_different_sig2 <-
function(model,method='Nelder-Mead',init.optim=NULL){
  if(is.null(init.optim)==T){init.optim=c(rep(-5,length(model$trees)),rep(0,model$ncoeff))}
  else{}
  opt=optim(par=init.optim,fn=model$fun,method=method,control=list(maxit=5000))
  par_fixed=model$par_fixed
  par_names=c() 
  for (i in 1:length(model$trees)){par_names=c(par_names,paste('sigsq_tree_',i,sep=''))}
  par_names=c(par_names,'a','b','c')
  par=list()
  for (i in 1:length(model$trees)){par[[i]]=2*exp(opt$par[i]) ; names(par)[i]=par_names[i]}
  if (model$ncoeff==3){par$a=opt$par[(length(model$trees)+1)] ; par$b=opt$par[(length(model$trees)+2)] ; par$c=opt$par[(length(model$trees)+3)]}
  else if (model$ncoeff==2){par$b=opt$par[(length(model$trees)+1)] ; par$c=opt$par[(length(model$trees)+2)]}
  else if (model$ncoeff==1){par$c=opt$par[(length(model$trees)+1)]}
  else{}
  # now retrieve x0: not needed actually, since this could be done using the ACE function
  return(res=list(lnL=-opt$value,aic=2*(length(init.optim)+opt$value),k=length(init.optim),par=par,par_fixed=par_fixed,convergence=opt$convergence,message=opt$message,trees=model$trees,traits=model$traits,Npts=model$Npts)) 
}
