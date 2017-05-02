Optim_bBM_x2x_flex_pts_multiple_starts <-
function(tree,trait,Npts=50,method='Nelder-Mead',verbose=T){
  pars=matrix(c(rep(log(var(trait)/(2*max(branching.times(tree)))),3),rep(1,3),rep(c(0,0.1,0.5),2)),nrow=6,ncol=2) # different starting points for sigsq and bounds
  starts=list()
  for (i in 1:dim(pars)[1]){
    starts[[i]]=try(Optim_bBM_x2x_flex_pts_start(tree, trait,Npts= Npts,method=method,start.point=c(pars[i,1],0,0,min(trait)-pars[i,2]*(max(trait)-min(trait)),max(trait)+ pars[i,2]*(max(trait)-min(trait)))))
    if (class(starts[[i]])=='try-error') {starts[[i]]=list(lnL=-Inf)}
    if (verbose==T){cat("Finished preliminary fit, log-lik= ",starts[[i]]$lnL,sep="\n")}
  }
  # add one fit with bounds fixed?
  starts[[(length(starts)+1)]]=try(Optim_bBM_quadratic(tree, trait,Npts= Npts,bounds=c(min(trait),max(trait)),method=method))
  if (class(starts[[length(starts)]])=='try-error') {starts[[length(starts)]]=list(lnL=-Inf)}
  if (verbose==T){cat("Finished preliminary fit with bounds fixed to min/max, log-lik= ",starts[[length(starts)]]$lnL,sep="\n")}
  lnls=c(starts[[1]]$lnL)
  for (i in 2:length(starts)){lnls=c(lnls,starts[[i]]$lnL)}
  mod=starts[[which(lnls==max(lnls))]] # the starting point which got us to the highest lnl
  # check for convergence before returning ML estimate
  if (verbose==T){cat("Starting final fit... ",sep="\n")}
  temp=try(Optim_bBM_x2x_flex_pts_start(tree, trait,Npts= Npts,method=method,start.point=c(log(mod$par$sigsq/2),mod$par$b,mod$par$c,mod$par$bounds[1],mod$par$bounds[2])))
  rerun=0
  if (class(temp)=='try-error'){
    temp=mod
    if (verbose==T){cat("Final optimization step produced an error: reporting the best of the initial optimization steps.")}
  }
  else {
    while((temp$convergence!=0)&(rerun<11)){
      if (verbose==T){cat("Trying to improve convergence... ",sep="\n")
        if (rerun==10){cat("...this is the last try!",sep="\n")}}
      temp=try(Optim_bBM_x2x_flex_pts_start(tree, trait,Npts= Npts,method= method,start.point=c(log(temp$par$sigsq/2),temp$par$b, temp$par$c, temp$par$bounds[1], temp$par$bounds[2])))
      rerun=rerun+1
      if (class(temp)=='try-error'){
        temp=mod ; rerun=20
        if (verbose==T){cat("Final optimization step produced an error: reporting the best of the initial optimization steps.")}
      }
    }
  }
  return(temp)
}
