LogLik_bounds <-
function(tree_formatted,dCoeff,dMat,bounds){
  #tree_formatted obtained through FormatTree ; dCoeff=log(sigsq/2)
  Npts=dim(dMat$diag)[1]
  tree_formatted2= tree_formatted
  pMat=prep_mat_exp(dCoeff,dMat,bounds) # edited
  logFactor=0
  for (i in 1:dim(tree_formatted2$tab)[1]){
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat = pMat)
    norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
    tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
    logFactor=logFactor+log(norm)
  }
  return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
  # this is where we take the max over the root position
}
