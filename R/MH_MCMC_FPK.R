MH_MCMC_FPK <-
function(tree,trait,bounds,Nsteps=500000,record_every=100,plot_every=500,Npts=50,pars_init=c(0,0,0,0,25),prob_update=c(0.2,0.2,0.2,0.2,0.2),verbose=TRUE,plot=TRUE,save_to='MCMC_FPK_test.Rdata',save_every=10000,type_priors=c(rep('Normal',4),'Uniform'),shape_priors=list(c(0,10),c(0,10),c(0,10),c(0,10),NA),proposal_type='Uniform',proposal_sensitivity=c(0.1,0.1,0.1,0.1,1),prior.only=F,burnin.plot=0.1){
# prior.only to sample from prior only (check that MCMC algorithm mixes well). Default to F for actual posterior exploration	
# burnin.plot gives the proportion burnin for plots only (the whole chain is actually saved)  
# we update parameters separately: prob_update gives the probability that each param is updated
  if (is.numeric(trait)){
    if ((min(trait)<bounds[1])|(max(trait)>bounds[2])){stop('Some values in the trait vector exceed the bounds.')} 
  }
  if (class(trait)=='list') {
    if ((min(unlist(trait))<bounds[1])|(max(unlist(trait))>bounds[2])){stop('Some values in the trait data exceed the bounds.')}
  }
SEQ=seq(from=-1.5,to=1.5,length.out= Npts) # the potential V is modelled as a quadratic function over [-1.5,1.5], but in real data space, this corresponds to [bounds[1],bounds[2]]
V_init= pars_init[2]*SEQ^4+pars_init[3]*SEQ^2+pars_init[4]*SEQ
temp= pars_init
chain=matrix(NA,Nsteps/record_every,11)
colnames(chain)=c('step','sigsq','a','b','c','root','lnprior','lnlik','quasi-lnpost','Accept','Par_updated')  
if (prior.only==T){lnlik=1}
else {
lnlik= LogLik_bounds_est_root(tree, trait,dCoeff=temp[1],x0_pos=temp[5],V= V_init,bounds=bounds)
}
lnprior= log_prior_5pars_root_bounds(type=type_priors,shape=shape_priors,pars=temp,Npts=Npts)
lnpost=lnlik+ lnprior
if ((is.na(lnpost))|(lnpost==(-Inf))){stop('Likelihood cannot be estimated at initial parameters. Please change them')}
for (i in 1:Nsteps){
	par_to_update=sample(1:length(pars_init),size=1,prob=prob_update) # sample which parameter will be updated in this step
	sensitivity_temp=rep(0,length(pars_init)) # set all sensitivities to 0 so that parameters are not updated...
	sensitivity_temp[par_to_update]= proposal_sensitivity[par_to_update] #... except the one chosen
	prop= proposal_5pars_root_bounds(type='Uniform',sensitivity= sensitivity_temp,pars=temp)
	lnprior_proposed= log_prior_5pars_root_bounds(type=type_priors,shape=shape_priors,pars=prop,Npts=Npts)
	if (lnprior_proposed ==(-Inf)){lnpost_proposed=-Inf} # no lnl calculation when prior is null
	else {
	if (prior.only==T){lnlik_proposed=1}
	else {
	#SEQ=seq(from=-1.5,to=1.5,length.out= Npts)		
	V_proposed= prop[2]*SEQ^4+ prop[3]*SEQ^2+ prop[4]*SEQ # proposed potential
	lnlik_proposed= try(LogLik_bounds_est_root(tree,trait, dCoeff=prop[1],x0_pos=prop[5], V=V_proposed,bounds=bounds)) # test with try
	}
	if (class(lnlik_proposed)=='try-error'){lnlik_proposed=NaN} # redo steps where likelihood calculation failed
	lnpost_proposed= lnlik_proposed + lnprior_proposed # un-normalized log-posterior
	}
	ALPHA=exp(lnpost_proposed-lnpost) # acceptance ratio (ratio of un-norm. posteriors)
	if (is.nan(ALPHA)){i=i-1} # re-do this step if it produces an NaN
	else {
	U=runif(1)
	if (U<ALPHA){
		temp= prop
		lnlik= lnlik_proposed
		lnpost= lnpost_proposed
		lnprior= lnprior_proposed
		accept=1
	}
	else {accept=0}
	if (i%%record_every==0){
	chain[(i/record_every),]=c(i,2*exp(temp[1]),temp[2:4],bounds[1]+(temp[5]-1)*(bounds[2]-bounds[1])/(Npts-1),lnprior,lnlik, lnpost,accept, par_to_update)
	}
	if (i%%plot_every==0){
	if (verbose==T){
		print(chain[(i/record_every),])
	}
	if (plot==T){
par(mfrow=c(3,3))
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),2],type='l',main='sigsq',log='y',ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),6],type='l',main='root',ylab='',xlab='')
if (is.numeric(trait)){
  abline(h=min(trait)+(max(trait)-min(trait))/2,col=2)
}
if (class(trait)=='list') {
  abline(h=min(unlist(trait))+(max(unlist(trait))-min(unlist(trait)))/2,col=2)
}
plot(1,1,ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),3],type='l',main='a (x^4 term)',ylab='',xlab='')
abline(h=0,col=2)
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),4],type='l',main='b (x^2 term)',ylab='',xlab='')
abline(h=0,col=2)
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),5],type='l',main='c (x term)',ylab='',xlab='')
abline(h=0,col=2)
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),7],type='l',main='lnprior',ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),8],type='l',main='lnlik',ylab='',xlab='')
plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),9],type='l',main='quasi-lnpost',ylab='',xlab='')
	}
	}
	if (i%%save_every==0){
		save(chain,file=save_to)
		}	
	}
}
return(chain)
}
