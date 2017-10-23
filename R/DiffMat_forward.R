DiffMat_forward <-
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
  eig=eigen(M)
  passage=matrix(NA,dim(eig$vectors)[1],dim(eig$vectors)[2])
  for (col in 1:dim(eig$vectors)[2]){passage[,col]=Re(eig$vectors[,col])}
  return(list(Diff=M,diag=diag(Re(eig$values)),passage=passage))
}
