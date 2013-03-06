#' \code{permutation} of the specified lmDME object
#'
#' Produces de specified lmDME plus the required permuted objects (sampling the
#' columns), using the same parameters to fit the additional models.
#' 
#' @param model formula object to carry out the decomposition.
#' @param data data.frame with individuals (rows) and samples/conditions (columns)
#' @param design data.frame with the design of the experiment, (rows) 
#'  samples/conditions as in data columns and as many columns to indicate the
#'  factors presents in each sample.
#' @param Bayes Should limma estimate empirical Bayes statistics, i. e., 
#'  moderated t-statistics? Default value is FALSE.
#' @param verbose Should the process progress be printed? Default value is FALSE.
#' @param NPermutations number of permutations to be calculated. Default value 
#'  is 100.
#' @param nCpus number of cores to be used. Default value is 1, i.e. sequential
#'  calculation.
#' @param ... Additional parameters for \code{\link{lmFit}} function.
#'
#' @return 
#'  \item{list}{contains the original lmDME object plus the required amount 
#'  of permuted versions.} 
#'
#' @seealso \code{\link{lmdme}}
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#'
#' @examples
#' {
#' data(stemHypoxia)
#' 
#' ##Just to make a balance dataset in the Fisher sense (2 samples per 
#' ## time*oxygen levels)
#' design<-design[design$time %in% c(0.5, 1, 5) & design$oxygen %in% c(1,5,21),]
#' design$time<-as.factor(design$time)
#' design$oxygen<-as.factor(design$oxygen)
#' rownames(M)<-M[, 1]
#' 
#' ##Keeping appropriate samples only
#' M<-M[, colnames(M) %in% design$samplename] 
#' 
#' ##Just to test if it works. On real scenario use NPermutations >= 100 if the
#' ## conditions (columns) of M allows it. Verbose parameter is FALSE by default
#' permuted<-permutation(model=~time*oxygen, data=M, design=design,
#' NPermutations=2, nCpus=3)
#' }
#'
#' @exportMethod permutation
#' @docType methods
#' @name permutation
#' @rdname lmDME-permutation
#' @aliases permutation-methods
setGeneric(name="permutation", def=function(model, data, design, Bayes=FALSE,
  verbose=FALSE, NPermutations=100, nCpus=1, ...){
  standardGeneric("permutation")
})
#'
#' @name permutation
#' @rdname lmDME-permutation
#' @inheritParams permutation
#' @usage \S4method{permutation}{formula,data.frame,data.frame}(model,data,design,Bayes=FALSE,verbose=FALSE,NPermutations=100,nCpus=1,...)
#' @aliases permutation,formula,data.frame,data.frame-method
setMethod(f="permutation", signature=signature(model="formula",
  data="data.frame", design="data.frame"), definition=function(model, data,
  design, Bayes=FALSE, verbose=FALSE, NPermutations=100, nCpus=1, ...){
  ##Generate the permuted samples: original structure + NPermutations
  permutations<-cbind(1:ncol(data),
    sapply(1:NPermutations, function(iteration, indexes){sample(indexes)},
    indexes=1:ncol(data)))
  
  ##Auxiliary functions, print if verbose==TRUE
  if(verbose){
    printnow<-function(...){ cat(...);flush.console()}
  }else{
    printnow<-function(...){invisible (NULL)}
  }

  ##Auxiliary functions for parallel processing if available
  if(require(parallel)){
    parlapply<-mclapply
    ##Get the cpus data for parallel lmdme calculation
    ##Check if windows platform
    if(.Platform$OS.type == "windows"){
      nCpus<-1
    }else{
      ncores<-detectCores()
      printnow("using", nCpus, "core/s from", ncores, "available/s\n")
    }
  }
  else{
    ##parallel not installed so, use the well known lapply
    parlapply<-function(X, FUN, ..., mc.cores){lapply(X, FUN, ...)}
  }

  ##Get the cpus data for parallel lmdme calculation
  ##Check if windows platform
  if(.Platform$OS.type == "windows"){
    nCpus<-1
  }else{
      ncores<-detectCores()
      printnow("using", nCpus, "core/s from", ncores, "available/s\n")
  }

  ##Just to make the index of the permutations
  permutedModels<-as.list(1:ncol(permutations))
  names(permutedModels)<-c("Original", as.character(1:NPermutations))

  ##Start calculating the permutations
  time <- Sys.time()
  printnow("Start Time:", as.character(time), "\n")
  permutedModels<-parlapply(permutedModels, function(index){
    printnow("Running Model: ", index, "\n")
    return(lmdme(model=model, data=data[, permutations[,index]], design=design,
      Bayes=Bayes, verbose=verbose, ...))}, mc.cores = nCpus)
  ##End calculations
  timeEnd <- Sys.time()
  printnow("End Time", as.character(timeEnd), "\n")
  timedifference<-timeEnd-time
  printnow("Time difference of ", format(unclass(timedifference),
    digits=getOption("digits")), " ", attr(timedifference , "units"), "\n")
  return(permutedModels)
})
