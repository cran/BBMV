proposal_nclades_plus_3_pars <-
function(type='Uniform',sensitivity,pars,n_clades){ 
  # same params: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c)
  npars=n_clades+3
  if (type=='Uniform'){
    par_temp=rep(NA,length(pars))
    for (i in 1:npars){
      par_temp[i]=runif(n=1,min=pars[i]-sensitivity[i],max=pars[i]+sensitivity[i])
    }
  }
  if (type=='Normal'){
    par_temp=rep(NA,length(pars))
    for (i in 1:npars){
      par_temp[i]=rnorm(n=1,mean=pars[i],sd= sensitivity[i])
    }
  }
  return(par_temp)
}
