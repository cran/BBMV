log_prior_7pars_root_bounds <-
function(type=c(rep('Normal',4),rep('Uniform',3)),shape=list(c(0,10),c(0,10),c(0,10),c(0,10),NA,10,10),pars,Npts_int,trait){
	# pars: the actual parameters for which to calculate the prior (dCoeff,a,b,c,x0,bmin,bmax). Although bmin and bmax only move between points of the grid, they are treated as continuous for practical purposes (all likelihood functions treat them as such...) 
	# type: either uniform of normal prior for each param, only uniform (discrete) for root position, and only uniform continuous for bounds
	# Npts_int is Npts in the interval, but if bounds lie further than min/max of the trait, then we have a real Npts that takes these extra points into account
	# the prior is on log(sigsq/2)=dCoeff, not sigsq
	# shape: params of the prior, min/max for uniform , mean/sd for normal, no shape for root position yet (uniform discrete between 1 and Npts)
	# shape[6,7]: the maximum number of steps explored away from the bounds
	step=(max(trait)-min(trait))/(Npts_int-1)
	Npts=floor((pars[7]-pars[6])/step+1) # to prevent rounding errors
	Npts_max= Npts_int+shape[[6]][1]+shape[[7]][1]
	p=list() 
	for (i in 1:4){
		if (type[i]=='Normal'){p[[i]]=dnorm(x=pars[i],mean=shape[[i]][1],sd=shape[[i]][2])}
		if (type[i]=='Uniform'){p[[i]]=dunif(x=pars[i],min=shape[[i]][1],max=shape[[i]][2])}
	}
	# root should also be treated as a continuous value???? change the call in log-lik function then to retrieve the exact point. Or we shift it when we update bmin... better probably
	if (type[5]=='Uniform'){
		if (pars[5]%in%c(1:Npts)){
			p[[5]]=1/Npts_max
		}
		else {p[[5]]=0}
	}
	# change that and make bounds only move on the grid I think
	if (type[6]=='Uniform'){p[[6]]=dunif(x=pars[6],min=min(trait)-shape[[6]][1]*step,max=min(trait))
		}
	if (type[7]=='Uniform'){p[[7]]=dunif(x=pars[7],min=max(trait),max=max(trait)+shape[[7]][1]*step)
		}	
	
# # 	else { stop('Only uniform prior on root position supported so far. More fancy things coming soon.')}
	return(log(p[[1]]*p[[2]]*p[[3]]*p[[4]]*p[[5]]*p[[6]]*p[[7]]))
}
