add.ML.landscape.FPK=function(fit,Npts=100,COLOR=1,LTY='dashed'){
  if ('a'%in%names(fit$par)){a=fit$par$a}
  else {a=fit$par_fixed$a}
  if ('b'%in%names(fit$par)){b=fit$par$b}
  else {b=fit$par_fixed$b}
  if ('c'%in%names(fit$par)){c=fit$par$c}
  else {c=fit$par_fixed$c}
  # build potential and stationary distribution of the trait
  bounds=fit$par_fixed$bounds
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  V=a*SEQ^4+b*SEQ^2+c*SEQ #potential
  step=(bounds[2]-bounds[1])/(Npts-1)
  lines((exp(-V)/sum(exp(-V)*step))~seq(from=bounds[1],to=bounds[2],length.out=Npts),col=COLOR,lty=LTY,lwd=3)
}
