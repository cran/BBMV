fit_BBMV <-
function(tree,trait,Npts=50,method='Nelder-Mead',verbose=T,V_shape,bounds='estimate',init.optim=NULL){
  if (bounds[1]=='estimate'){
    if ((V_shape%in%c('flat','linear','quadratic','full')==F)){
      stop("Wrong specification of V_shape: should be one of 'flat','linear','quadratic', or 'full'")
    }
    if (V_shape=='flat'){
      return(Optim_bBM_0_flex_pts_multiple_starts(tree,trait,Npts,method=method,verbose=verbose))
    }
    if (V_shape=='linear'){
      return(Optim_bBM_x_flex_pts_multiple_starts(tree,trait,Npts,method=method,verbose=verbose))
    }
    if (V_shape=='quadratic'){
      return(Optim_bBM_x2x_flex_pts_multiple_starts(tree,trait,Npts,method=method,verbose=verbose))
    }
    if (V_shape=='full'){
      return(Optim_bBM_x4x2x_flex_pts_multiple_starts(tree,trait,Npts,method=method,verbose=verbose))
    }
  }
  if ((class(bounds)=='numeric')&(length(bounds)==2)){
    if (class(V_shape)=="character"){
      if ((V_shape%in%c('flat','linear','quadratic','full')==F)){
        stop("Wrong specification of V_shape: should be one of 'flat','linear','quadratic', or 'full'")
      }
      if (V_shape=='flat'){
        return(Optim_bBM_bounds_fixed_potential(tree,trait,V=rep(0,Npts),bounds=bounds))
      }
      if (V_shape=='linear'){
        return(Optim_bBM_linear(tree,trait,Npts,bounds=bounds,method=method,init.optim=init.optim))
      }
      if (V_shape=='quadratic'){
        return(Optim_bBM_quadratic(tree,trait,Npts,bounds=bounds,method=method,init.optim=init.optim))
      }
      if (V_shape=='full'){
        return(Optim_bBM_x4x2x(tree,trait,Npts,bounds=bounds,method=method,init.optim=init.optim))
      }
    }
    if (class(V_shape)=="numeric"){
      print("Ignoring argument 'Npts' and using the length of the 'V_shape' vector instead.")
      return(Optim_bBM_bounds_fixed_potential(tree,trait,V=V_shape,bounds=bounds,init.optim=init.optim))
    }
    else {stop("Argument 'V_shape' should be either a character or a numeric vector of values for the potential")}
  }
  else {stop("Argument 'bounds' should be either a numeric vector specifying the two bounds that are fixed or be set to 'estimate'")
  }
}
