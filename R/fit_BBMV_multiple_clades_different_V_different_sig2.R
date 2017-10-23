fit_BBMV_multiple_clades_different_V_different_sig2 <-
function(trees,traits,bounds,a=NULL,b=NULL,c=NULL,Npts=50,method='Nelder-Mead',init.optim=NULL){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  lnls=ks=rep(NA,length(trees))
  fits=list()
  for (i in 1:length(trees)){
    lnl_temp=lnL_BBMV(trees[[i]],traits[[i]],bounds=bounds,a=a,b=b,c=c,Npts=Npts)
    fit_temp=find.mle_FPK(model=lnl_temp,method=method,init.optim=init.optim)
    fits[[i]]=fit_temp ; names(fits)[i]=paste('fit_clade_',i,sep='')
    lnls[i]=fit_temp$lnL ; ks[i]=fit_temp$k
  }
  return(list(lnL=sum(lnls),aic=2*(sum(ks)-sum(lnls)),k=sum(ks),fits=fits))
}
