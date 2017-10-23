log_prior_nclades_plus_3_pars <-
function(type=NULL,shape=NULL,pars,n_clades){
  # params: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c)
  # type: either uniform of normal prior for each param, only uniform (discrete) for root position
  # the prior is on log(sigsq/2)=dCoeff, not sigsq
  # shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
  npars=n_clades+3
  if (is.null(type)){
    type=list()
    shape=list()
    for (i in 1:npars){type[[i]]='Normal'; shape[[i]]=c(0,10)}
  }
  p=rep(NA,npars)
  for (i in 1:npars){
    if (type[i]=='Normal'){p[i]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
    if (type[i]=='Uniform'){p[i]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
  }
  return(sum(log(p)))
}
