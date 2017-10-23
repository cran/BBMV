MH_MCMC_FPK_multiclades <-
function(trees,traits,bounds,Nsteps=500000,record_every=100,plot_every=500,Npts=50,pars_init=NULL,prob_update=NULL,verbose=TRUE,plot=TRUE,save_to='MCMC_FPK_test.Rdata',save_every=10000,type_priors=NULL,shape_priors=NULL,proposal_type='Normal',proposal_sensitivity=NULL,prior.only=F,burnin.plot=0.1){
  # the oder of parameters is the same for pars_init, prob_update,type_priors,shape_priors and proposal_sensitivity. It is: (dCoeff1,dCoeff2,...,dCoeffn,a,b,c) , with dCoeff=log(sigsq/2)
  # prior.only to sample from prior only (check that MCMC algorithm mixes well). Default to F for actual posterior exploration	
  # burnin.plot gives the proportion burnin for plots only (the whole chain is actually saved)  
  # we update parameters separately: prob_update gives the probability that each param is updated
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  n_clades=length(trees) ; npars=n_clades+3
  for (i in 1:n_clades){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
    ###### new piece of code added for Measurment error incorporation
    if (is.numeric(traits[[i]])){
      if ((min(traits[[i]])<bounds[1])|(max(traits[[i]])>bounds[2])){stop(paste('Some values in trait ',i,' vector exceed the bounds.'))} 
    }
    if (class(traits[[i]])=='list') {
      if ((min(unlist(traits[[i]]))<bounds[1])|(max(unlist(traits[[i]]))>bounds[2])){stop(paste('Some values in trait ',i,' vector exceed the bounds.'))}
    }
  }  
  if (is.null(prob_update)){prob_update=rep(1/npars,npars)}
  if (is.null(pars_init)){pars_init=c(rnorm(n=n_clades,mean=-8,sd=3),rnorm(n=3,mean=0,sd=2))}
  if (is.null(proposal_sensitivity)){proposal_sensitivity=rep(0.1,npars)}
  ######  end new code  
  SEQ=seq(from=-1.5,to=1.5,length.out= Npts) # the potential V is modelled as a quadratic function over [-1.5,1.5], but in real data space, this corresponds to [bounds[1],bounds[2]]
  V_init= pars_init[(n_clades+1)]*SEQ^4+pars_init[(n_clades+2)]*SEQ^2+pars_init[(n_clades+3)]*SEQ
  temp= pars_init
  chain=matrix(NA,Nsteps/record_every,(n_clades+9))
  colnames(chain)[c(1,(n_clades+2):(n_clades+9))]=c('step','a','b','c','lnprior','lnlik','quasi-lnpost','Accept','Par_updated') 
  for (clade in 1:n_clades){eval(parse(text=paste("colnames(chain)[",clade,"+1]='sigsq_clade_",clade,"'",sep='')))}
  if (prior.only==T){lnlik=1}
  else {
    NEG_LNL_func=lnl_BBMV_multiclades_same_V_different_sig2(trees=trees,traits=traits,bounds=bounds,a=NULL,b=NULL,c=NULL,Npts=50)$fun
    LNL_func=function(X){return(-NEG_LNL_func(X))}
    lnlik= LNL_func(X=temp)
  }
  lnprior= log_prior_nclades_plus_3_pars(type=type_priors,shape=shape_priors,pars=temp,n_clades=n_clades)
  lnpost=lnlik+ lnprior
  if ((is.na(lnpost))|(lnpost==(-Inf))){stop('Likelihood cannot be estimated at initial parameters. Please change them')}
  for (i in 1:Nsteps){
    par_to_update=sample(1:length(pars_init),size=1,prob=prob_update) # sample which parameter will be updated in this step
    sensitivity_temp=rep(0,length(pars_init)) # set all sensitivities to 0 so that parameters are not updated...
    sensitivity_temp[par_to_update]= proposal_sensitivity[par_to_update] #... except the one chosen
    prop= proposal_nclades_plus_3_pars(type=proposal_type,sensitivity=sensitivity_temp,pars=temp,n_clades=n_clades)
    lnprior_proposed=log_prior_nclades_plus_3_pars(type=type_priors,shape=shape_priors,pars=prop,n_clades=n_clades) 
    if (lnprior_proposed ==(-Inf)){lnpost_proposed=-Inf} # no lnl calculation when prior is null
    else {
      if (prior.only==T){lnlik_proposed=1}
      else {
        #	V_proposed= prop[(n_clades+1)]*SEQ^4+ prop[(n_clades+2)]*SEQ^2+ prop[(n_clades+3)]*SEQ # proposed potential
        lnlik_proposed= try(LNL_func(X=prop))
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
        chain[(i/record_every),c(1,(n_clades+2):(n_clades+9))]=c(i,temp[(n_clades+1):(n_clades+3)],lnprior,lnlik, lnpost,accept, par_to_update)
        for (clade in 1:n_clades){chain[(i/record_every),(clade+1)]=2*exp(temp[clade])}
      }
      if (i%%plot_every==0){
        if (verbose==T){
          print(chain[(i/record_every),])
        }
        if (plot==T){
          par(mfrow=c(ceiling(sqrt(n_clades+6)),ceiling(sqrt(n_clades+6))))
          for (clade in 1:n_clades){plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(clade+1)],type='l',main=paste('sigsq_clade_',clade,sep=''),log='y',ylab='',xlab='')}
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+2)],type='l',main='a (x^4 term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+3)],type='l',main='b (x^2 term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+4)],type='l',main='c (x term)',ylab='',xlab='')
          abline(h=0,col=2)
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+5)],type='l',main='lnprior',ylab='',xlab='')
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+6)],type='l',main='lnlik',ylab='',xlab='')
          plot(chain[floor((i/record_every)*burnin.plot):(i/record_every),(n_clades+7)],type='l',main='quasi-lnpost',ylab='',xlab='')
        }
      }
      if (i%%save_every==0){
        save(chain,file=save_to)
      }	
    }
  }
  return(chain)
}
