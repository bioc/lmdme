#' Getters for lmDME object
#'
#' Obtain lmDME slot information, acording to the given function call (see Values). If term parameter is not specified, it will return all the available terms, otherwise just the one specified.
#' 
#' @param object lmDME class object.
#' @param term character with the corresponding term/s to return. Default value is NULL in order to return every available term/s.
#' 
#' @return according to the call one of the following objects can be returned
#'  \item{design}{used experiment design data.frame.}
#'  \item{model}{used decompose formula.}
#'  \item{modelDecomposition}{list of decomposed model formulas.}
#'  \item{residuals or resid, coef or coefficients, fitted or fitted.values, p.values or F.p.values}{list of appropieate slot where each item is a matrix that will have G rows(individuals) x k columns(levels of the corresponding model term).}
#'  \item{components}{list with corresponding PCA or PLSR term according to decomposition function call.}
#'  \item{componentsType}{character name vector with the information of the components calculation.}
#'  
#' @seealso \code{\link{lmdme}}, \code{\link{decomposition}}, \code{\link{print}}, \code{\link{show}}
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
#' fit.model  <- model(fit) #The model formula used
#' fit.design <- design(fit)#The design data.frame used
#' fit.modelDecomposition<- modelDecomposition(fit) #How the decomposition was carried out
#' ##Getting the coefficients
#' timeCoef <- coef(fit,term="time")  #or coef(fit) to get all term coefficients 
#' fit.p.values <- p.values(fit,term="time") #for the same term
#' fit.f.values <- F.p.values(fit,term="time") #or the F-values 
#' ##Getting the residuals or fitted values
#' interactionResid <- resid(fit, term="time:oxygen") #or resid(fit) to get all term residuals
#' oxygenDFit <- fitted(fit, term="(Intercept)") #or fitted(fit) to get all term fitted values
#' 
#' @exportMethod fitted.values
#' @docType methods
#' @name fitted.values
#' @rdname lmDME-getters
#' @importFrom stats fitted.values
#' @usage \S4method{fitted.values}{lmDME}(object,term)
#' @aliases fitted.values,lmDME-method
setMethod(f="fitted.values",signature = "lmDME", definition = function(object, term=NULL){  
    ##If Term == NULL the full decomposed fitted.values
    if(is.null(term)){term<-names(object@coefficients)}
    ##Check for term in model
    if(!any(term %in% names(object@coefficients))){stop("ERROR: ", toString(term), " not in model: ", object@model)}
    
    ##Get the fitted values for the requested term/s
    output <- lapply(as.list(term),function(termino){
      ##The intercept is a special case, for the rest is the term-1 for the model matrix
      ##Then hat(y) = X(the model matrix) %*% hat(beta)
      if(termino == "(Intercept)"){
	mm <- model.matrix(as.formula(paste("~ 1",sep="")),data=object@design)
      }
      else{
	mm <- model.matrix(as.formula(paste("~",termino,"-1",sep="")),data=object@design)
      }
      return(t(apply(object@coefficients[[termino]],MARGIN=1,FUN=function(x,mm){mm%*%x},mm)))
    })
    
    names(output)<-term
    return(output)
})
#'
#' @exportMethod fitted
#' @name fitted
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats fitted
#' @usage \S4method{fitted}{lmDME}(object,term)
#' @aliases fitted,lmDME-method
setMethod(f="fitted",signature = "lmDME", definition = function(object, term=NULL){return(fitted.values(object,term))})
#'
#' @exportMethod coef
#' @name coef
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats coef
#' @usage \S4method{coef}{lmDME}(object,term)
#' @aliases coef,lmDME-method
setMethod(f="coef",signature = "lmDME", definition = function(object, term=NULL){
  ##If Term == NULL the full decomposed coefficients
  if(is.null(term)){return(object@coefficients)}
  ##Else search for the term over the decomposed coefficients
  if(all(term %in% names(object@coefficients))){return(object@coefficients[term])}else{stop("ERROR: ", toString(term), " not in model: ", object@model)}
})
#'
#' @exportMethod coefficients
#' @name coefficients
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats coefficients
#' @usage \S4method{coefficients}{lmDME}(object,term)
#' @aliases coefficients,lmDME-method
setMethod(f="coefficients",signature = "lmDME", definition = function(object, term=NULL){return(coef(object,term))})
#'
#' @exportMethod resid
#' @name resid
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats resid
#' @usage \S4method{resid}{lmDME}(object,term)
#' @aliases resid,lmDME-method
setMethod(f="resid",signature = "lmDME", definition = function(object, term=NULL){
  ##If Term == NULL the full decomposed residuals
  if(is.null(term)) return(object@residuals)
  ##Else search for the term over the decomposed residuals  
  if(all(term %in% names(object@residuals))){return(object@residuals[term])}else{stop("ERROR: ", toString(term), " not in model: ", object@model)}
})
#'
#' @exportMethod residuals
#' @name residuals
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @importFrom stats residuals
#' @usage \S4method{residuals}{lmDME}(object,term)
#' @aliases residuals,lmDME-method
setMethod(f="residuals",signature = "lmDME", definition = function(object, term=NULL){return(resid(object,term))})
#'
#' @exportMethod F.p.values
#' @name F.p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases F.p.values
setGeneric(name="F.p.values",def = function(object, term=NULL){standardGeneric("F.p.values")})
#'
#' @name F.p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases F.p.values,lmDME-method
setMethod(f="F.p.values",signature = "lmDME", definition = function(object, term=NULL){
  ##If Term == NULL the full decomposed F.p.values    
  if(is.null(term)) return(object@F.p.values)
  ##Else search for the term over the F.p.values  
  if(all(term %in% names(object@F.p.values))){return(object@F.p.values[term])}else{stop("ERROR: ", toString(term), " not in model: ", object@model)}
})
#'
#' @exportMethod p.values
#' @name p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases p.values,lmDME-method
setGeneric(name="p.values",def = function(object, term=NULL){standardGeneric("p.values")})
#'
#' @name p.values
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases p.values,lmDME-method
setMethod(f="p.values",signature = "lmDME", definition = function(object, term=NULL){
  ##If Term == NULL the full decomposed p.values
  if(is.null(term)) return(object@p.values)
  ##Else search for the term over the p.values  
  if(all(term %in% names(object@p.values))){return(object@p.values[term])}else{stop("ERROR: ", toString(term), " not in model: ", object@model)}
})
#'
#' @exportMethod modelDecomposition
#' @name modelDecomposition
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases modelDecomposition,lmDME-method
setGeneric(name="modelDecomposition",def = function(object, term=NULL){standardGeneric("modelDecomposition")})
#'
#' @name modelDecomposition
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases modelDecomposition,lmDME-method
setMethod(f="modelDecomposition",signature = "lmDME", definition = function(object, term=NULL){
  ##If Term == NULL the full decomposed models
  if(is.null(term)){return(object@modelDecomposition)}
  ##Else search for the term over the p.values  
  if(all(term %in% names(object@modelDecomposition))){return(object@modelDecomposition[term])}else{stop("ERROR: ", toString(term), " not in model: ", object@model)}
})
#'
#' @exportMethod components
#' @name components
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases components,lmDME-method
setGeneric(name="components",def = function(object, term=NULL){standardGeneric("components")})
#'
#' @name components
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases components,lmDME-method
setMethod(f="components",signature = "lmDME", definition = function(object, term=NULL){
  ##If Term == NULL the full decomposed models
  if(is.null(term)){return(object@components)}
  ##Else search for the term over the p.values  
  if(all(term %in% names(object@components))){return(object@components[term])}else{stop("ERROR: ", toString(term), " not in model/s: ", toString(names(object@components)))}
})
#'
#' @exportMethod componentsType
#' @name componentsType
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases componentsType,lmDME-method
setGeneric(name="componentsType",def = function(object, term=NULL){standardGeneric("componentsType")})
#'
#' @name componentsType
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases componentsType,lmDME-method
setMethod(f="componentsType",signature = "lmDME", definition = function(object, term=NULL){return(object@componentsType)})
#'
#' @exportMethod model
#' @name model
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases model,lmDME-method
setGeneric(name="model",def = function(object){standardGeneric("model")})
#'
#' @name model
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases model,lmDME-method
setMethod(f="model",signature = "lmDME", definition = function(object){return(object@model)})
#'
#' @exportMethod design
#' @name design
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases design,lmDME-method
setGeneric(name="design",def = function(object){standardGeneric("design")})
#'
#' @name design
#' @rdname lmDME-getters
#' @inheritParams fitted.values
#' @aliases design,lmDME-method
setMethod(f="design",signature = "lmDME", definition = function(object){return(object@design)})

