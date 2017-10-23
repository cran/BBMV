fit_FPK_multiple_clades_different_V_different_sig2 <-
function(trees,traits,a=NULL,b=NULL,c=NULL,Npts=50,method='Nelder-Mead',init.optim=NULL){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_BBMV instead')}
  # define bounds far way  
  # min_trait=min(traits[[1]]) ; max_trait=max(traits[[1]])
  # for (i in 2:length(trees)){
  #   min_trait=min(min_trait,min(traits[[i]])) ; max_trait=max(max_trait,max(traits[[i]]))
  # }
  # bounds=c(min_trait-(max_trait-min_trait)/2,max_trait+(max_trait-min_trait)/2)
  # new formulation for when there is measurement error
  span=c(min(unlist(traits)),max(unlist(traits)))
  bounds=c(span[1]-(span[2]-span[1])/2,span[2]+(span[2]-span[1])/2)  
  return(fit_BBMV_multiple_clades_different_V_different_sig2(trees=trees,traits=traits,bounds=bounds,a=a,b=b,c=c,Npts=Npts,method=method,init.optim=init.optim))
}
