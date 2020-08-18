parallel_rf <- function(mod_formula, data, ntree=500, n.cores=NA){
  #load libraries.
  library(randomForest)
  library(doParallel)
  
  #register parallel environment.
  if(is.na(n.cores)){
    n.cores <- detectCores()
  }
  registerDoParallel(cores = n.cores)
  
  #fit model.
  to_return <- foreach(i = 1:n.cores) %dopar% {
    return(randomForest(mod_formula, data = data, ntree=round(ntree/n.cores)))
  }
  to_return <- do.call(combine, to_return)
  
  #return output.
  return(to_return)
}