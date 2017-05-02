FormatTree_bounds <-
function(tree,trait,V,bounds){
  tree=reorder.phylo(tree,'postorder')
  ntips=length(tree$tip.label)
  tab=cbind(tree$edge,tree$edge.length) ; colnames(tab)=c('parent','children','brlen')
  Pos=list() # one element per node
  for (i in 1:(2*ntips-1)){
    if (i>ntips){
      Pos[[i]]=1
    }
    else {
      Pos[[i]]= VectorPos_bounds(trait[tree$tip.label[i]],V,bounds=bounds)
    }
  }
  return(list(tab=tab,Pos=Pos))
}
