VectorPos_bounds <-
function(x,V,bounds){
  Npts=length(V)
  if (length(x)==1){ # only one value per tip
    X=rep(0,Npts)  
    if (x==bounds[2]){X[Npts]=1}
    else {
      nx=(Npts-1)*(x-bounds[1])/(bounds[2]-bounds[1])
      ix=floor(nx)
      ux=nx-ix
      X[ix+2]=ux
      X[ix+1]=1-ux
    }	
  }
  else {
    # here we treat the case in which we do not have a single value but a vector of values measured
    MAT=matrix(0,length(x),Npts)
    for (i in 1:length(x)){ # problem with the recursive form here
      if (x[i]==bounds[2]){MAT[i,Npts]=1}
      else {
        nx=(Npts-1)*(x[i]-bounds[1])/(bounds[2]-bounds[1])
        ix=floor(nx)
        ux=nx-ix
        MAT[i,ix+2]=ux
        MAT[i,ix+1]=1-ux
      }	
    }
    X=apply(MAT,2,mean) # and average them
  }
  return(X*(Npts-1)/(bounds[2]-bounds[1]))
}
