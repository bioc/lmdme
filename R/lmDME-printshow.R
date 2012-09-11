#' \code{Show}, \code{Print}  or \code{Summary} a lmDME object
#'
#' Generic Show/Print/Summary Method for lmDME class output visualization.
#'
#' @param x lmDME class object
#' @param object lmDME class object
#' @param term character with the corresponding term to return. Default value is NULL to return every decomposed term (if more than one available)
#'
#' @return according to the call
#'  \item{print, show or summary:}{console output text with increasing detail of lmDME object.}
#'  \item{show or summary}{console output text of the lmDME object, plus a data.frame with model decomposition summary data.}
#'
#' @seealso \code{\link{lmdme}}, \code{\link{coef}}, \code{\link{resid}}, \code{\link{fitted}}, \code{\link{modelDecomposition}}, \code{\link{components}}, \code{\link{componentsType}}
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#'
#' @examples
#' data(stemHypoxia)
#' ##Just to make a balance dataset in the Fisher sense (2 samples per time*oxygen levels)
#' design<-design[design$time %in% c("0.5","1","5") & design$oxygen %in% c("1","5","21"),]
#' design$time <-as.factor(design$time)
#' design$oxygen<-as.factor(design$oxygen)
#' rownames(M)<-M[,1]
#' M <- M[,colnames(M) %in% design$samplename] #Keeping appropiate samples only
#' ##ANOVA decomposition
#' fit <- lmdme(model=~time+oxygen+time:oxygen,data=M,design=design)
#' \dontrun{ 
#'   fit #equivalent to call show(fit) 
#'   print(fit)
#'   summary(fit) 
#' }
#'
#' @exportMethod print
#' @docType methods
#' @name print
#' @rdname lmDME-printshow
#' @usage \S4method{print}{lmDME}(x,term)
#' @aliases print,lmDME-method
setMethod(f="print",signature="lmDME", definition = function(x, term=NULL){
   show(x)
   ##If term is available the data of the corresponding term
   if(!missing(term)){
    cat("Residuals (head):\n"); print(head(residuals(x,term=term)))
    cat("Coefficients (head):\n"); print(head(coefficients(x,term=term)))
    cat("P-Values (head):\n"); print(head(pvalues(x,term=term)))
  }else{
    ##All the terms of the model
    cat("Residuals (head):\n"); print(lapply(residuals(x,term=term),head))
    cat("Coefficients (head):\n"); print(lapply(coefficients(x,term=term),head))
    cat("P-Values (head):\n"); print(lapply(pvalues(x,term=term),head))
  }
})
#' @exportMethod show
#' @name show
#' @rdname lmDME-printshow
#' @inheritParams print
#' @usage \S4method{show}{lmDME}(object)
#' @aliases show,lmDME-method
setMethod(f="show",signature="lmDME", definition = function(object){
    ##Basic information of lmDME class
    cat("lmDME object: \n")
    cat("Data dimension: ", nrow(object@residuals[[1]])," x ", ncol(object@residuals[[1]]))
    cat("\nDesign (head): \n"); print(head(object@design))
    cat("\nModel:");print(object@model)
    cat("Model deflaction: \n"); 
    deflaction<-data.frame(Step = 1:length(object@modelDecomposition), Names = names(object@modelDecomposition),Formula = c("~ 1", paste("~ -1 + ",names(object@modelDecomposition)[-1],sep="")),CoefCols = as.numeric(lapply(object@coefficients,ncol)))
    print(deflaction)
   
    ##Just to keep de information of how the deflaction process was carried out
    return(invisible(deflaction))
})
#'
#' @exportMethod summary 
#' @name summary
#' @rdname lmDME-printshow
#' @inheritParams show
#' @usage \S4method{summary}{lmDME}(object)
#' @aliases summary,lmDME-method
setMethod(f="summary",signature = "lmDME", definition = function(object){return(show(object))})


