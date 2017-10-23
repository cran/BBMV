get.landscape.FPK.MCMC <-
function(chain,bounds,Npts=100,burnin=0.1,probs.CI=c(0.05,0.95),COLOR_MEDIAN='red',COLOR_FILL='red',transparency=0.3,main='Macroevolutionary landscapes MCMC',ylab='N.exp(-V)',xlab='Trait',xlim=NULL,ylim=NULL){
  chain2=chain[-c(1:floor(dim(chain)[1]*burnin)),]  # remove burnin
  all_V=as.data.frame(matrix(NA,dim(chain2)[1],Npts))
  step=(bounds[2]-bounds[1])/(Npts-1)
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  for (i in 1:dim(chain2)[1]){
    temp=chain2[i,'a']*SEQ^4+chain2[i,'b']*SEQ^2+chain2[i,'c']*SEQ #potential
    all_V[i,]=exp(-temp)/sum(exp(-temp)*step)
  }
  if (is.null(xlim)){xlim=bounds}
  if (is.null(ylim)){ylim=c(0,max(all_V))}
  par(mfrow=c(1,1))
  plot(apply(all_V,2,function(x){quantile(x,probs=c(0.5))})~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=COLOR_MEDIAN,lwd=3,type='l',main=main,ylab=ylab,xlab=xlab,xlim=xlim,ylim=ylim)
  polygon(x=c(seq(from=bounds[1],to=bounds[2],length.out=100),seq(from=bounds[2],to=bounds[1],length.out=100)),y=c(apply(all_V,2,function(x){quantile(x,probs=probs.CI[1])}),rev(apply(all_V,2,function(x){quantile(x,probs=probs.CI[2])}))),col=adjustcolor(col=COLOR_FILL,alpha.f=transparency))
}
