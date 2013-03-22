#' \code{leverage} test of lmdme objects
#' 
#' This function calculates the leverage test for each individual using the 
#' Principal Component Analysis (comps function) on the coefficients of the
#' given decomposed model term.
#'
#' @param object lmdme class object.
#' @param comps a numeric vector indicating the PCA component indexes to keep. 
#'  Default the first two components (1:2).
#' @param term a character specifying the model term. 
#' @param level the quantile level. Default value 0.95
#'
#' @return data.frame with the following fields
#'  \item{leverage}{numeric for the corresponding row leverage}
#'  \item{over}{logical indicating if the leverage > quantile(leverage,level) 
#'  for the given decomposed term}
#'
#' @seealso \code{\link{prcomp}}, \code{\link{quantile}}
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#'
#' @references Tarazona S, Prado-Lopez S, Dopazo J, Ferrer A, Conesa A, 
#'  Variable Selection for Multifactorial Genomic Data, Chemometrics and
#'  Intelligent Laboratory Systems, 110:113-122 (2012)  
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
#' ##Leverages for the first two Principal Components and q95 (default value). 
#' ##Leverages for the first three Principal Components and q99.
#' leverages2PCDefault<-leverage(fit, term="time:oxygen")
#' leverages3PCq99<-leverage(fit, comps=1:3, term="time:oxygen", level=0.99)
#' }
#'
#' @exportMethod leverage
#' @docType methods
#' @name leverage
#' @rdname lmdme-leverage
#' @aliases leverage-methods
setGeneric(name="leverage", def=function(object, comps=1:2, term, level=0.95){
  standardGeneric("leverage")
})
#'
#' @name leverage
#' @rdname lmdme-leverage
#' @inheritParams leverage
#' @aliases leverage,lmdme-method
setMethod(f="leverage", signature="lmdme", definition=function(object,
  comps=1:2, term, level=0.95){
  ##Check for term presence
  if(missing(term)){
    stop("ERROR: missing term")
  }
  ##Obtain the leverages
  pca<-prcomp(coefficients(object, term=term))
  out<-data.frame(leverage=apply(pca$x[, comps]^2, 1, sum))
  out$over<-out$leverage > quantile(out$leverage, level)
  return(out)
})
