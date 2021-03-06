lnL_FPK <-
function(tree,trait,a=NULL,b=NULL,c=NULL,Npts){
  if (is.numeric(trait)){
    bounds=c(min(trait)-(max(trait)-min(trait))/2,max(trait)+(max(trait)-min(trait))/2)}
  if (class(trait)=='list') {
    bounds=c(min(unlist(trait))-(max(unlist(trait))-min(unlist(trait)))/2,max(unlist(trait))+(max(unlist(trait))-min(unlist(trait)))/2)
  }
  if (sum(tree$tip.label%in%names(trait))<max(length(trait),length(tree$tip.label))){stop('Tip names in tree do not match names of the trait vector')}
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  tree_formatted= FormatTree_bounds(tree,trait,V=rep(0,Npts),bounds=bounds)
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T)
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0)  
      if (is.null(c)==F){fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds))}}
      # only c varies (e.g. flat landscape if a=b=0)  
      else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[2]*SEQ),bounds=bounds))}}
    }
    # a is fixed (e.g. quadratic landscape if a=0)  
    else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(a*SEQ^4+X[2]*SEQ^2+X[3]*SEQ),bounds=bounds))}}
  }
  # the full model: no parameter fixed
  else {fun=function(X){return(-LogLik_bounds(tree_formatted=tree_formatted,dCoeff=X[1],dMat=DiffMat_backwards(X[2]*SEQ^4+X[3]*SEQ^2+X[4]*SEQ),bounds=bounds))}}
  return(list(fun=fun,ncoeff=ncoeff,par_fixed=list(a=a,b=b,c=c,bounds=bounds),tree=tree,trait=trait,Npts=Npts))
}
