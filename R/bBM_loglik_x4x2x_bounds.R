bBM_loglik_x4x2x_bounds <-
function(tree_formatted,Npts=100,bounds){ 
  fun=function(X){ # X=c(dCoeff,a,b,c) No intercept since we are only interested in V'
    SEQ=seq(from=-1.5,to=1.5,length.out=Npts) # very sensitive
    V=X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ
    dMat=DiffMat_backwards(V)
    return(-LogLik_bounds(tree_formatted,X[1],dMat,bounds))
  }
  return(fun)
}
