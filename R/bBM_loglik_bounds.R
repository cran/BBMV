bBM_loglik_bounds <-
function(tree_formatted,dMat,bounds){
  fun=function(dCoeff){
    return(-LogLik_bounds(tree_formatted,dCoeff,dMat,bounds))
  }
  return(fun)
}
