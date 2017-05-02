bBM_loglik_quadra_bounds <-
function(tree_formatted,Npts=100,bounds){ 
  fun=function(X){ # X=c(dCoeff,c,b) where c is the linear term for the potential, b is the quadratic. No intercept since we are only interested in V'
    SEQ=seq(from=0,to=1,length.out=Npts)
    V=X[2]*SEQ^2+X[3]*SEQ # corrected: order of coeeficients matches function with estimated bounds
    dMat=DiffMat_backwards(V)
    return(-LogLik_bounds(tree_formatted,X[1],dMat,bounds))
  }
  return(fun)
}
