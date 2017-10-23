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
      if (class(trait)=='numeric'){Pos[[i]]= VectorPos_bounds(trait[tree$tip.label[i]],V=V,bounds=bounds)} # only one value per tip
      else {Pos[[i]]= VectorPos_bounds(trait[[tree$tip.label[i]]],V=V,bounds=bounds)} # multiple values per tip (i.e. uncertainty)
    }
  }
  return(list(tab=tab,Pos=Pos))
}
