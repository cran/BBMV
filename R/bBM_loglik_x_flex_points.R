bBM_loglik_x_flex_points <-
function(tree,trait,Npts){ 
  fun=function(X){ # X=c(dCoeff,c,bmin,bmax) where c is the linear term for the potential. No intercept since we are only interested in V'
    # Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
    step=(max(trait)-min(trait))/(Npts-1)
    steps_below=ceiling((min(trait)-X[3])/step) 
    bmin=min(trait)-steps_below*step # we approximate X[3] by the nearest point on the grid, closer to min(trait)
    steps_above=floor((X[4]-max(trait))/step)
    bmax=max(trait)+steps_above*step # we approximate X[4] by the nearest point on the grid, closer to max(trait)
    Npts_tot=Npts+ steps_below+ steps_above
    SEQ=seq(from=-1.5,to=1.5,length.out=Npts_tot) # very sensitive
    V=X[2]*SEQ
    return(-LogLik_bounds_est(tree,trait,X[1],V,c(X[3],X[4])))
  }
  return(fun)
}
