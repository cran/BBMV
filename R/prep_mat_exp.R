prep_mat_exp <-
function(dCoeff,dMat,bounds){
  vDiag=dMat$diag ; P=dMat$passage ; tP=solve(P,tol = 1e-30) #; tP=t(P)
  Npts=dim(dMat$diag)[1]
  tau=((bounds[2]-bounds[1])/(Npts-1))^2
  diag_expD=exp(dCoeff)/tau*diag(vDiag) # faster than the for loop
  return(list(P=P,tP=tP,diag_expD=diag_expD))
}
