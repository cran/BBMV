bBM_loglik_0_flex_points <-
function(tree,trait,Npts){ 
  fun=function(X){ # X=c(dCoeff,bmin,bmax) 
    # Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
    step=(max(trait)-min(trait))/(Npts-1)
    steps_below=ceiling((min(trait)-X[2])/step) 
    bmin=min(trait)-steps_below*step # we approximate X[2] by the nearest point on the grid, closer to min(trait)
    steps_above=floor((X[3]-max(trait))/step)
    bmax=max(trait)+steps_above*step # we approximate X[3] by the nearest point on the grid, closer to max(trait)
    Npts_tot=Npts+ steps_below+ steps_above
    V=rep(0,Npts_tot)
    return(-LogLik_bounds_est(tree,trait,X[1],V,c(X[2],X[3])))
  }
  return(fun)
}
