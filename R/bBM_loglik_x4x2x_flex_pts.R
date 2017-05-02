bBM_loglik_x4x2x_flex_pts <-
function(tree,trait,Npts){ 
  fun=function(X){ # X=c(dCoeff,a(x4),b(x2),c(x),bmin,bmax) No intercept since we are only interested in V'
    # Npts is the number of points IN THE TRAIT INTERVAL: WE ADD POINTS OUTSIDE WHEN FURTHERING THE BOUNDS
    step=(max(trait)-min(trait))/(Npts-1)
    steps_below=ceiling((min(trait)-X[5])/step) 
    bmin=min(trait)-steps_below*step # we approximate X[5] by the nearest point on the grid, closer to min(trait)
    steps_above=floor((X[6]-max(trait))/step)
    bmax=max(trait)+steps_above*step # we approximate X[6] by the nearest point on the grid, closer to max(trait)
    Npts_tot=Npts+ steps_below+ steps_above
    SEQ=seq(from=-1.5,to=1.5,length.out=Npts_tot) # very sensitive
    V=X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ
    return(-LogLik_bounds_est(tree,trait,X[1],V,c(bmin,bmax)))
  }
  return(fun)
}
