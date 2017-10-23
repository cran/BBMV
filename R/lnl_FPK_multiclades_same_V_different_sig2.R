lnl_FPK_multiclades_same_V_different_sig2 <-
function(trees,traits,a=NULL,b=NULL,c=NULL,Npts=50){
  if (length(trees)!=length(traits)){stop('The list of trees and the list of traits differ in length.')}
  if (length(trees)==1){stop('There is only one tree and trait vector: use the function lnl_FPK instead')}
  for (i in 1:length(trees)){
    if (sum(trees[[i]]$tip.label%in%names(traits[[i]]))<max(length(traits[[i]]),length(trees[[i]]$tip.label))){stop(paste('Tip names in tree ',i,' do not match names of corresponding trait vector'))}
  }
# define bounds far way  
#   min_trait=min(traits[[1]]) ; max_trait=max(traits[[1]])
# for (i in 2:length(trees)){
#   min_trait=min(min_trait,min(traits[[i]])) ; max_trait=max(max_trait,max(traits[[i]]))
# }
#   bounds=c(min_trait-(max_trait-min_trait)/2,max_trait+(max_trait-min_trait)/2)
# new formulation for when there is measurement error
  span=c(min(unlist(traits)),max(unlist(traits)))
  bounds=c(span[1]-(span[2]-span[1])/2,span[2]+(span[2]-span[1])/2)
  SEQ=seq(from=-1.5,to=1.5,length.out=Npts)
  trees_formatted=list()
  for (i in 1:length(trees)){
    trees_formatted[[i]]=FormatTree_bounds(trees[[i]],traits[[i]],V=rep(0,Npts),bounds=bounds)
  }
  ncoeff=(is.null(a)==T)+(is.null(b)==T)+(is.null(c)==T) # OK
  npar=ncoeff+length(trees)
  par_names=c() 
  for (i in 1:length(trees)){par_names=c(par_names,paste('dCoeff_tree_',i,sep=''))}
  par_names=c(par_names,'a','b','c')
  if (is.null(a)==F){
    if (is.null(b)==F){
      # all three shape parameters fixed (e.g. flat landscape if a=b=c=0): seems to work 
      if (is.null(c)==F){
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+c*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }      
      
      # only c varies (e.g. flat landscape if a=b=0): seems to work   
      else {
        fun_text='fun=function(X){return('
        for (i in 1:length(trees)){
          fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+b*SEQ^2+X[',length(trees)+1,']*SEQ),bounds=bounds)',sep='') 
        }
        fun_text=paste(fun_text,')}',sep='') # the end parenthesis
        fun=eval(parse(text=fun_text))
      }
    }
    # only a is fixed (e.g. quadratic landscape if a=0): seems to work 
    else {
      fun_text='fun=function(X){return('
      for (i in 1:length(trees)){
        fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(a*SEQ^4+X[',length(trees)+1,']*SEQ^2+X[',length(trees)+2,']*SEQ),bounds=bounds)',sep='') 
      }
      fun_text=paste(fun_text,')}',sep='') # the end parenthesis
      fun=eval(parse(text=fun_text))
    } 
  }
  
  # the full model: no parameter fixed
  else {
    fun_text='fun=function(X){return('
    for (i in 1:length(trees)){
      fun_text=paste(fun_text,'-LogLik_bounds(tree_formatted=trees_formatted[[',i,']],dCoeff=X[',i,'],dMat=DiffMat_backwards(X[',length(trees)+1,']*SEQ^4+X[',length(trees)+2,']*SEQ^2+X[',length(trees)+3,']*SEQ),bounds=bounds)',sep='') 
    }
    fun_text=paste(fun_text,')}',sep='') # the end parenthesis
    fun=eval(parse(text=fun_text))
  }
  return(list(fun=fun,ncoeff=ncoeff,npar=npar,par_names=par_names,par_fixed=list(a=a,b=b,c=c,bounds=bounds),trees=trees,traits=traits,Npts=Npts))
}
