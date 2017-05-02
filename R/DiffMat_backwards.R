DiffMat_backwards <-
function (V){
  # V is a vector representing the potential, with 'Npts' numeric values
  Npts=length(V)
  M=matrix(0,Npts,Npts)
  for (i in 2:(Npts-1)){
    M[i-1,i]=exp((V[i]-V[i-1])/2)
    M[i+1,i]=exp((V[i]-V[i+1])/2)
    M[i,i]=-(M[i-1,i]+M[i+1,i])
  }
  M[2,1]=exp((V[1]-V[2])/2)
  M[1,1]=-M[2,1]
  M[Npts-1,Npts]=exp((V[Npts]-V[Npts-1])/2)
  M[Npts,Npts]=-M[Npts-1,Npts]
  M=t(M) # we go backwards in time!
  eig=eigen(M)
  return(list(Diff=M,diag=diag(eig$values),passage=eig$vectors))
}
