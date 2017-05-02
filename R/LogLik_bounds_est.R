LogLik_bounds_est <-
function(tree,trait ,dCoeff,V,bounds){
  if ((bounds[1]>min(trait))|(bounds[2]<max(trait))) {return(-Inf)} # bounds have to be outside the trait interval
  else {
    Npts_tot=length(V)		
    tree_formatted=FormatTree_bounds(tree,trait,V,bounds)		
    dMat=DiffMat_backwards(V)
    pMat=prep_mat_exp(dCoeff,dMat,bounds) # edited
    tree_formatted2= tree_formatted
    logFactor=0
    
    for (i in 1:dim(tree_formatted2$tab)[1]){
      tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]*ConvProp_bounds(X=tree_formatted2$Pos[[tree_formatted2$tab[i,2]]],t=tree_formatted2$tab[i,3],prep_mat=pMat)
      norm=sum(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]])
      tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]= tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]/norm
      logFactor=logFactor+log(norm)
    }
    return(log(max(tree_formatted2$Pos[[tree_formatted2$tab[i,1]]]))+logFactor)
    # this is where we take the max over the root position
  }
}
