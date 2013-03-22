#' \code{p.adjust} of p-values for Multiple Test Comparison Corrections
#'
#' Given a set of p-values, returns adjusted p-values using one of several
#' methods.
#' 
#' @param p numeric vector of p-values as in stats::p.adjust or lmdme class
#' object.
#' @param method correction method available in \code{\link{p.adjust.methods}}.
#' @param term character with the corresponding term to return.
#' @param ... other arguments.
#' @param drop should try to drop the list structure if length==1? Default value
#' is TRUE 
#'
#' @return according to the call, one of the following objects can be returned
#' \item{numeric}{vector of adjusted p-values.}
#' \item{matrix}{for lmdme object If term!=NULL, the corresponding character is
#'  looked up within the list of p.values returning the associated matrix of G
#'  rows (individuals) x k columns (levels of the corresponding model term) with
#'  the adjusted p-values.}
#'  
#' @seealso \code{\link{p.adjust}}, \code{\link{p.adjust.methods}}
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#' 
#' @examples
#' {
#' data(stemHypoxia)
#' 
#' ##Just to make a balanced dataset in the Fisher sense (2 samples per 
#' ## time*oxygen levels)
#' design<-design[design$time %in% c(0.5, 1, 5) & design$oxygen %in% c(1,5,21),]
#' design$time<-as.factor(design$time)
#' design$oxygen<-as.factor(design$oxygen)
#' rownames(M)<-M[, 1]
#' 
#' ##Keeping appropriate samples only
#' M<-M[, colnames(M) %in% design$samplename] 
#' 
#' ##ANOVA decomposition
#' fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)
#' 
#' ##Adjust p-values only on the interaction p.values using false discovery rate
#' ## method
#' pInteraction<-p.values(fit, term="time:oxygen") 
#' FDRInteraction<-p.adjust(fit, term="time:oxygen", method="fdr")
#' corrected<-sum(pInteraction < 0.05) - sum(FDRInteraction < 0.05)
#' }
#'
#' @exportMethod p.adjust
#' @docType methods
#' @name p.adjust
#' @rdname lmdme-padjust
#' @aliases p.adjust-methods
setGeneric(name="p.adjust", def=function(p, ...){standardGeneric("p.adjust")})
#'
#' @exportMethod p.adjust
#' @name p.adjust
#' @rdname lmdme-padjust
#' @inheritParams p.adjust
#' @aliases p.adjust,ANY-method
setMethod(f="p.adjust", signature="ANY", definition=stats::p.adjust)
#'
#' @exportMethod p.adjust
#' @name lmdme-padjust
#' @rdname lmdme-padjust
#' @inheritParams p.adjust
#' @usage \S4method{p.adjust}{lmdme}(p, term=NULL, method=p.adjust.methods, drop=TRUE)
#' @aliases p.adjust,lmdme-method
setMethod(f="p.adjust", signature="lmdme", definition=function(p, term=NULL,
  method=p.adjust.methods, drop=TRUE){
  out<-lapply(p.values(p, term, drop=FALSE), function(x){
    apply(x, MARGIN=2, FUN=p.adjust, method=method)
  })

  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
