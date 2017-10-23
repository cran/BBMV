reformat_multiclade_results <-
function(fit){
  n_clades=length(fit$trees)
  fits=list()
  if ('sigsq_tree_1' %in% names(fit$par)){ # then we have different sigmas in each clade
    for (i in 1:n_clades){
      eval(parse(text=paste('fits$fit_clade_',i,'=list()',sep='')))
      fits[[i]]$par=list()
      fits[[i]]$par$sigsq=fit$par[[i]]
      if ('a'%in%names(fit$par)){fits[[i]]$par$a=fit$par$a}
      if ('b'%in%names(fit$par)){fits[[i]]$par$b=fit$par$b}
      if ('c'%in%names(fit$par)){fits[[i]]$par$c=fit$par$c}
      fits[[i]]$par_fixed=fit$par_fixed
      fits[[i]]$tree=fit$trees[[i]]
      fits[[i]]$trait=fit$trait[[i]]
      fits[[i]]$Npts=fit$Npts
    }
  } 
  else {
    for (i in 1:n_clades){
      eval(parse(text=paste('fits$fit_clade_',i,'=list()',sep='')))
      fits[[i]]$par=list()
      fits[[i]]$par$sigsq=fit$par$sigsq
      if ('a'%in%names(fit$par)){fits[[i]]$par$a=fit$par$a}
      if ('b'%in%names(fit$par)){fits[[i]]$par$b=fit$par$b}
      if ('c'%in%names(fit$par)){fits[[i]]$par$c=fit$par$c}
      fits[[i]]$par_fixed=fit$par_fixed
      fits[[i]]$tree=fit$trees[[i]]
      fits[[i]]$trait=fit$trait[[i]]
      fits[[i]]$Npts=fit$Npts
    }  
  }
  return(fits)
}
