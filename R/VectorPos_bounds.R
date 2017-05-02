VectorPos_bounds <-
function(x,V,bounds){
  Npts=length(V)
  X=rep(0,Npts)
  if (x==bounds[2]){X[Npts]=1}
  else {
    nx=(Npts-1)*(x-bounds[1])/(bounds[2]-bounds[1])
    ix=floor(nx)
    ux=nx-ix
    X[ix+2]=ux
    X[ix+1]=1-ux
  }	
  return(X*(Npts-1)/(bounds[2]-bounds[1]))
}
