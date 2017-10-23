trans_to_fixed <-
function(x,bounds){
  return(1.5*(2*(x-bounds[1])/(bounds[2]-bounds[1])-1))
}
