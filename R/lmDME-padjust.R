#' \code{p.adjust} of p-values for Multiple Test Comparisons Correction
#'
#' Given a set of p-values, returns p-values adjusted using one of several methods.
#' 
#' @param p numeric vector of p-values as in stats::p.adjust or lmDME class object.
#' @param method correction method availables in p.adjust.methods.
#' @param term character with the corresponding term to return.
#' @param ... other arguments.
#'
#' @return according to the call one of the following objects can be returned
#' \item{numeric}{vector of adjusted p-values}
#' \item{matrix}{for lmDME object If term!=NULL, the corresponding character is looked up within list of p.values returning the associated matrix of G rows(individuals) x k columns(levels of the corresponding model term) with the adjusted p-values.}
#'  
#' @seealso \code{\link{p.adjust}}, \code{\link{p.adjust.methods}}
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
#' ##p.adjust only over interction p.values using false descovery rate method
#' pInteraction <- p.values(fit,term="time:oxygen")[[1]] 
#' FDRIntercept <- p.adjust(fit,term="time:oxygen",method="fdr")[[1]]
#' corrected <- sum(pInteraction < 0.05) - sum(FDRIntercept < 0.05)
#' 
#' @exportMethod p.adjust
#' @docType methods
#' @name p.adjust
#' @rdname lmDME-padjust
#' @aliases p.adjust-methods
setGeneric(name="p.adjust",def = function(p, ...){standardGeneric("p.adjust")})
#'
#' @exportMethod p.adjust
#' @name p.adjust
#' @rdname lmDME-padjust
#' @inheritParams p.adjust
#' @aliases p.adjust,ANY-method
setMethod(f="p.adjust",signature = "ANY", definition = stats::p.adjust)
#'
#' @exportMethod p.adjust
#' @name lmDME-padjust
#' @rdname lmDME-padjust
#' @inheritParams p.adjust
#' @usage \S4method{p.adjust}{lmDME}(p,term,method)
#' @aliases p.adjust,lmDME-method
setMethod(f="p.adjust",signature = "lmDME", definition = function(p,term=NULL,method = p.adjust.methods)
{ return(lapply(p.values(p,term),function(x){apply(x,MARGIN=2,FUN=p.adjust,method=method)}))
})