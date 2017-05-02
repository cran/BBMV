bBM_loglik_linear_bounds <-
function(tree_formatted,Npts=100,bounds){ 
  fun=function(X){ # X=c(dCoeff,a) where a is the linear term for the potential. No intercept since we are only interested in V'
    SEQ=seq(from=0,to=1,length.out=Npts)
    V=X[2]*SEQ
    dMat=DiffMat_backwards(V)
    return(-LogLik_bounds(tree_formatted,X[1],dMat,bounds))
  }
  return(fun)
}
