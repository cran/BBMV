ConvProp_bounds <-
function(X,t,prep_mat){
  Npts=length(X)
  expD=matrix(0,Npts,Npts)
  diag(expD)=exp(t*prep_mat$diag_expD)
  a=prep_mat$P%*%expD%*%prep_mat$tP%*%X
  return(apply(a,1,function(x) max(x,0))) # prevent rounding errors for small numbers
}
