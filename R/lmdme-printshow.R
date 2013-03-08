#' \code{Show}, \code{Print} or \code{Summary} a lmdme object
#'
#' Generic Show/Print/Summary method for lmdme class output visualization.
#'
#' @param x lmdme class object.
#' @param object lmdme class object.
#' @param term character with the corresponding term to return. Default value 
#'  is NULL to return every decomposed term (if more than one available).
#'
#' @return according to the call
#'  \item{print, show or summary:}{console output text with increasing detail of
#'  lmdme object.}
#'  \item{show or summary}{console output text of the lmdme object, plus a 
#'  data.frame with model decomposition summary data.}
#'
#' @seealso \code{\link{lmdme}}, \code{\link{coef}}, \code{\link{resid}}, 
#' \code{\link{fitted}}, \code{\link{modelDecomposition}},
#' \code{\link{components}}, \code{\link{componentsType}}
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
#' ##ANOVA decomposition
#' fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)
#' }
#' \dontrun{ 
#'  #equivalent to call show(fit) 
#'  fit 
#'  print(fit)
#'  summary(fit) 
#' }
#'
#' @exportMethod print
#' @docType methods
#' @name print
#' @rdname lmdme-printshow
#' @usage \S4method{print}{lmdme}(x, term)
#' @aliases print,lmdme-method
setMethod(f="print", signature="lmdme", definition=function(x, term=NULL){
   show(x)
   ##If term is available the data of the corresponding term
   if(!missing(term)){
    cat("Residuals (head):\n")
    print(head(residuals(x, term=term)))
    cat("Coefficients (head):\n")
    print(head(coefficients(x, term=term)))
    cat("P-Values (head):\n")
    print(head(pvalues(x, term=term)))
  }else{
    ##All the terms of the model
    cat("Residuals (head):\n")
    print(lapply(residuals(x, term=term, drop=FALSE), head))
    cat("Coefficients (head):\n")
    print(lapply(coefficients(x, term=term, drop=FALSE), head))
    cat("P-Values (head):\n")
    print(lapply(pvalues(x, term=term, drop=FALSE), head))
  }
})
#' @exportMethod show
#' @name show
#' @rdname lmdme-printshow
#' @inheritParams print
#' @usage \S4method{show}{lmdme}(object)
#' @aliases show,lmdme-method
setMethod(f="show", signature="lmdme", definition=function(object){
    ##Basic information of lmdme class
    cat("lmdme object: \n")
    cat("Data dimension: ", nrow(object@residuals[[1]])," x ", 
      ncol(object@residuals[[1]]))
    cat("\nDesign (head): \n")
    print(head(object@design))
    cat("\nModel:")
    print(object@model)
    cat("Model decomposition: \n")
    deflaction<-data.frame(Step=1:length(object@modelDecomposition), 
      Names=names(object@modelDecomposition), Formula=c("~ 1", 
      paste("~ -1 + ", names(object@modelDecomposition)[-1], sep="")),
      CoefCols=as.numeric(lapply(object@coefficients, ncol)))
    print(deflaction)
   
    ##Just to keep de information of how the deflaction process was carried out
    return(invisible(deflaction))
})
#'
#' @exportMethod summary 
#' @name summary
#' @rdname lmdme-printshow
#' @inheritParams show
#' @usage \S4method{summary}{lmdme}(object)
#' @aliases summary,lmdme-method
setMethod(f="summary", signature="lmdme", definition=function(object){
  return(show(object))
})
