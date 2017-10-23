find.mle_FPK_multiple_clades_same_V_same_sig2 <-
function(model,method='Nelder-Mead',init.optim=NULL){
    if(is.null(init.optim)==T){init.optim=c(-5,rep(0,model$ncoeff))}
    else{}
    opt=optim(par=init.optim,fn=model$fun,method=method,control=list(maxit=50000))
  par_fixed=model$par_fixed
  par=list()
  par$sigsq=2*exp(opt$par[1])
  if (model$ncoeff==3){par$a=opt$par[2] ; par$b=opt$par[3] ; par$c=opt$par[4]}
  else if (model$ncoeff==2){par$b=opt$par[2] ; par$c=opt$par[3]}
  else if (model$ncoeff==1){par$c=opt$par[2]}
  else{}
  # now retrieve x0: not needed actually, since this could be done using the uncertainty function
  return(res=list(lnL=-opt$value,aic=2*(length(init.optim)+opt$value),k=length(init.optim),par=par,par_fixed=par_fixed,convergence=opt$convergence,message=opt$message,trees=model$trees,traits=model$traits,Npts=model$Npts)) 
}
