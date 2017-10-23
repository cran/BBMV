log_prior_5pars_root_bounds <-
function(type=c(rep('Normal',4),'Uniform'),shape=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),pars,Npts){
	# pars: the actual parameters for which to calculate the prior (dCoeff,a,b,c,x0).
	# type: either uniform of normal prior for each param, only uniform (discrete) for root position
	# the prior is on log(sigsq/2)=dCoeff, not sigsq
	# shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
	p=list() 
	for (i in 1:4){
		if (type[i]=='Normal'){p[[i]]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
		if (type[i]=='Uniform'){p[[i]]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
	}
	# root should also be treated as a continuous value???? change the call in log-lik function then to retrieve the exact point. Or we shift it when we update bmin... better probably
	if (type[5]=='Uniform'){
		if (pars[5]%in%c(1:Npts)){
			p[[5]]=1/Npts
		}
		else {p[[5]]=0}
	}
	
 	else { stop('Only uniform prior on root position is supported.')}
	return(log(p[[1]]*p[[2]]*p[[3]]*p[[4]]*p[[5]]))
}
