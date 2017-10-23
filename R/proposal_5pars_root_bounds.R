proposal_5pars_root_bounds <-
function(type='Uniform',sensitivity,pars){ 
	# same 5 params: (dCoeff,a,b,c,x0)
	# for the root we draw a point on the grid nearby from root-sensitivity, to root+sensitivity, excluding root (i.e. we HAVE to move). Probably better to move one step at the time, i.e. sensitivity[5]=1. Same for the bounds
	if (type=='Uniform'){
		par_temp=c(runif(n=1,min=pars[1]-sensitivity[1],max=pars[1]+sensitivity[1]),runif(n=1,min=pars[2]-sensitivity[2],max=pars[2]+sensitivity[2]),runif(n=1,min=pars[3]-sensitivity[3],max=pars[3]+sensitivity[3]),runif(n=1,min=pars[4]-sensitivity[4],max=pars[4]+sensitivity[4]))
	}
	if (type=='Normal'){
	par_temp=c(rnorm(n=1,mean=pars[1],sd= sensitivity[1]),rnorm(n=1,mean=pars[2],sd= sensitivity[2]),rnorm(n=1,mean=pars[3],sd= sensitivity[3]),rnorm(n=1,mean=pars[4],sd= sensitivity[4]))
	}
	# update the root
	if (sensitivity[5]==0){root=pars[5]}
	else{
		root=pars[5]+sample(size=1,x=c(seq(from=-sensitivity[5],to=-1,by=1),seq(from=1,to=sensitivity[5],by=1))) }
	return(c(par_temp,root))
}
