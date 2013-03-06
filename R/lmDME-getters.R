#' Getters for lmDME object
#'
#' Obtain lmDME slot information, according to the given function call (see
#' Values). If term parameter is not specified, it will return all the available
#' terms. Otherwise, just the one specified.
#' 
#' @param object lmDME class object.
#' @param term character with the corresponding term/s to return. Default value
#' is NULL in order to return every available term/s.
#' @param drop should try to drop list structure if length==1? Default value
#' is TRUE 
#' 
#' @return according to the call one of the following objects can be returned
#'  \item{design}{used experiment design data.frame.}
#'  \item{model}{used decompose formula.}
#'  \item{modelDecomposition}{list of decomposed model formulas.}
#'  \item{residuals, resid, coef, coefficients, fitted, fitted.values, 
#'  p.values or F.p.values}{list of appropriate slot where each item is a matrix
#'  that will have G rows (individuals) x k columns(levels of the corresponding
#'  model term).}
#'  \item{components}{list with corresponding PCA or PLSR term according to 
#'  decomposition function call.}
#'  \item{componentsType}{character name vector with the information of the 
#'  components calculation.}
#'  
#' @seealso \code{\link{lmdme}}, \code{\link{decomposition}}, 
#'  \code{\link{print}}, \code{\link{show}}
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
#' 
#' ##Let's inspect how the decomposition process was carried out:
#' ##a) The model formula used
#' ##b) The design data.frame used 
#' ##c) The decomposition itself
#' fit.model<-model(fit) 
#' fit.design<-design(fit)
#' fit.modelDecomposition<-modelDecomposition(fit)
#' 
#' ##Getting the specific "time" term coefficients, p-values or F-values.
#' ## Omit "term" parameter to get available data for all terms.
#' timeCoef<-coef(fit,term="time")  
#' fit.p.values<-p.values(fit,term="time") 
#' fit.f.values<-F.p.values(fit,term="time") 
#' 
#' ##Getting the residuals or fitted values, for the interaction "time:oxygen"
#' ## term. Omit "term" parameter to get available data for all terms.
#' interactionResid<-resid(fit, term="time:oxygen") 
#' interactionFit<-fitted(fit, term="time:oxygen")
#' }
#' 
#' @exportMethod fitted.values
#' @docType methods
#' @name fitted.values
#' @rdname lmDME-getters
#' @importFrom stats fitted.values
#' @usage \S4method{fitted.values}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases fitted.values,lmDME-method
setMethod(f="fitted.values", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){  
  ##If Term == NULL the full decomposed fitted.values
  if(is.null(term)){
    term<-names(object@coefficients)
  }else{
    ##Check for term in model
    if(!any(term %in% names(object@coefficients))){
      stop("ERROR: ", toString(term), " not in model: ", object@model)
    }
  }
    
  ##Get the fitted values for the requested term/s
  output<-lapply(as.list(term), function(termino){
    ##The intercept is a special case, for the rest is the term-1 for the model
    ##matrix. Then hat(y) = X(the model matrix) %*% hat(beta)
    if(termino == "(Intercept)"){
      mm<-model.matrix(as.formula(paste("~ 1",sep="")), data=object@design)
    }else{
      mm<-model.matrix(as.formula(paste("~", termino, "-1", sep="")),
        data=object@design)
    }
      return(t(apply(object@coefficients[[termino]], MARGIN=1, FUN=function(x,
        mm){mm%*%x}, mm)))
    })
    
  names(output)<-term

  ##Check for drop parameter
  if(drop & length(output)==1){
    output<-output[[1]]
  }
  
  return(output)
})
#'
#' @exportMethod fitted
#' @name fitted
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats fitted
#' @usage \S4method{fitted}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases fitted,lmDME-method
setMethod(f="fitted", signature = "lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  return(fitted.values(object, term, drop))
})
#'
#' @exportMethod coef
#' @name coef
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats coef
#' @usage \S4method{coef}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases coef,lmDME-method
setMethod(f="coef", signature = "lmDME", definition = function(object,
  term=NULL, drop=TRUE){
  ##If Term == NULL the full decomposed coefficients
  if(is.null(term)){
    out<-object@coefficients
  }else{
    ##Else search for the term over the decomposed coefficients
    if(all(term %in% names(object@coefficients))){
      out<-object@coefficients[term]
    }else{
      stop("ERROR: ", toString(term), " not in model: ", object@model)
    }
  }

  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
#'
#' @exportMethod coefficients
#' @name coefficients
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats coefficients
#' @usage \S4method{coefficients}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases coefficients,lmDME-method
setMethod(f="coefficients", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  return(coef(object, term, drop))
})
#'
#' @exportMethod resid
#' @name resid
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats resid
#' @usage \S4method{resid}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases resid,lmDME-method
setMethod(f="resid", signature = "lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  ##If Term == NULL the full decomposed residuals
  if(is.null(term)){
    out<-object@residuals
  }else{
    ##Else search for the term over the decomposed residuals  
    if(all(term %in% names(object@residuals))){
      out<-object@residuals[term]
    }else{
      stop("ERROR: ", toString(term), " not in model: ", object@model)
    }
  }

  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
#'
#' @exportMethod residuals
#' @name residuals
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats residuals
#' @usage \S4method{residuals}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases residuals,lmDME-method
setMethod(f="residuals", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  return(resid(object, term, drop))
})
#'
#' @exportMethod F.p.values
#' @name F.p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage F.p.values(object, term=NULL, drop=TRUE)
#' @aliases F.p.values-methods
setGeneric(name="F.p.values", def=function(object, term=NULL, drop=TRUE){
  standardGeneric("F.p.values")
})
#'
#' @name F.p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{F.p.values}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases F.p.values,lmDME-method
setMethod(f="F.p.values", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  ##If Term == NULL the full decomposed F.p.values    
  if(is.null(term)){
    out<-object@F.p.values
  }else{
    ##Else search for the term over the F.p.values  
    if(all(term %in% names(object@F.p.values))){
      out<-object@F.p.values[term]
    }else{
      stop("ERROR: ", toString(term), " not in model: ", object@model)
    }
  }
  
  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
#'
#' @exportMethod p.values
#' @name p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage p.values(object, term=NULL, drop=TRUE)
#' @aliases p.values-methods
setGeneric(name="p.values", def=function(object, term=NULL, drop=TRUE){
  standardGeneric("p.values")
})
#'
#' @name p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{p.values}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases p.values,lmDME-method
setMethod(f="p.values", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  ##If Term == NULL the full decomposed p.values
  if(is.null(term)){
    out<-object@p.values
  }else{
    ##Else search for the term over the p.values  
    if(all(term %in% names(object@p.values))){
      out<-object@p.values[term]
    }else{
      stop("ERROR: ", toString(term), " not in model: ", object@model)
    }
  }

  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
#'
#' @exportMethod modelDecomposition
#' @name modelDecomposition
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage modelDecomposition(object, term=NULL, drop=TRUE)
#' @aliases modelDecomposition-methods
setGeneric(name="modelDecomposition", def=function(object, term=NULL,
  drop=TRUE){
  standardGeneric("modelDecomposition")
})
#'
#' @name modelDecomposition
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{modelDecomposition}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases modelDecomposition,lmDME-method
setMethod(f="modelDecomposition", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  ##If Term == NULL the full decomposed models
  if(is.null(term)){
    out<-object@modelDecomposition
  }else{
    ##Else search for the term over the p.values  
    if(all(term %in% names(object@modelDecomposition))){
      out<-object@modelDecomposition[term]
    }else{
      stop("ERROR: ", toString(term), " not in model: ", object@model)
    }
  }

  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
#'
#' @exportMethod components
#' @name components
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage components(object, term=NULL, drop=TRUE)
#' @aliases components-methods
setGeneric(name="components", def=function(object, term=NULL, drop=TRUE){
  standardGeneric("components")
})
#'
#' @name components
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{components}{lmDME}(object, term=NULL, drop=TRUE)
#' @aliases components,lmDME-method
setMethod(f="components", signature="lmDME", definition=function(object,
  term=NULL, drop=TRUE){
  ##If Term == NULL the full decomposed models
  if(is.null(term)){
    out<-object@components
  }else{
    ##Else search for the term over the p.values  
    if(all(term %in% names(object@components))){
      out<-object@components[term]
    }else{
      stop("ERROR: ", toString(term), " not in model/s: ",
        toString(names(object@components)))
    }
  }
  
  ##Check for drop parameter
  if(drop & length(out)==1){
    out<-out[[1]]
  }

  return(out)
})
#'
#' @exportMethod componentsType
#' @name componentsType
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage componentsType(object)
#' @aliases componentsType-methods
setGeneric(name="componentsType", def=function(object){
  standardGeneric("componentsType")
})
#'
#' @name componentsType
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{componentsType}{lmDME}(object)
#' @aliases componentsType,lmDME-method
setMethod(f="componentsType", signature="lmDME", definition=function(object){
  return(object@componentsType)
})
#'
#' @exportMethod model
#' @name model
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage model(object)
#' @aliases model-methods
setGeneric(name="model", def=function(object){standardGeneric("model")})
#'
#' @name model
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{model}{lmDME}(object)
#' @aliases model,lmDME-method
setMethod(f="model", signature="lmDME", definition=function(object){
  return(object@model)
})
#'
#' @exportMethod design
#' @name design
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage design(object)
#' @aliases design-methods
setGeneric(name="design", def=function(object){standardGeneric("design")})
#'
#' @name design
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @usage \S4method{design}{lmDME}(object)
#' @aliases design,lmDME-method
setMethod(f="design", signature="lmDME", definition=function(object){
  return(object@design)
})

