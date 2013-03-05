#' \code{loadingplot} of interaction PCA decomposed lmDME object
#' 
#' This function plots the PCA loadings for a given interaction (A:B) lmDME
#' object's components slot, for the given "pc" component. The user can choose
#' which term (A or B) is used for x-axis and y-axis functions (B or A)
#' respectively.
#'
#' @param object lmDME class object.
#' @param term.x,term.y character indicating the model principal factor for the
#'  interaction term (term.x:term.y or term.y:term.x) for the corresponding x or
#'  y axis.
#' @param pc integer indicating which principal component loading is to be 
#'  plotted on the y-axis. Default value is 1.
#' @param col which color to use for each level present in term.y. 
#' @param ord.x numeric indicating the term.x levels order, for plotting
#'  purposes. If missing the levels order is used.
#' @param ... additional parameters for matplot.
#'
#' @return loading plot of the selected interaction (term.x:term.y) components
#'  lmDME object's slot, if PCA decomposition was applied.
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#'
#' @examples
#' {
#' data(stemHypoxia)
#' 
#' ##Just to make a balance dataset in the Fisher sense (2 samples per 
#' ## time*oxygen levels) 
#' design<-design[design$time %in% c(0.5,1,5) & design$oxygen %in% c(1,5,21), ]
#' design$time<-as.factor(design$time)
#' design$oxygen<-as.factor(design$oxygen)
#' rownames(M)<-M[, 1]
#' 
#' #Keeping appropriate samples only
#' M<-M[, colnames(M) %in% design$samplename] 
#' 
#' ##ANOVA decomposition
#' fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)
#' 
#' ##ASCA for all the available terms, over those subjects/genes where at least
#' ##one interaction coefficient is statistically different from zero (F-test
#' ##over the coefficients).
#' id<-F.p.values(fit, term="time:oxygen")[[1]]<0.001
#' decomposition(fit, decomposition="pca", scale="row", subset=id) 
#'
#' \dontrun{
#'  loadingplot(fit, term.x="time", term.y="oxygen") 
#'
#'  ##Or change the axis order
#'  loadingplot(fit, term.x="oxygen", term.y="time")
#'
#'  ##Or change the PC to display
#'  loadingplot(fit, term.x="time", term.y="oxygen", pc=2)
#'
#'  ##Or the order of x-levels
#'  loadingplot(fit, term.x="time", term.y="oxygen", ord.x=3:1)
#' }
#' }
#'
#' @exportMethod loadingplot
#' @importFrom pls loadingplot 
#' @docType methods
#' @name loadingplot
#' @rdname lmDME-loadingplot
#' @usage \S4method{loadingplot}{lmDME}(object, term.x, term.y, pc=1, ord.x, col, ...)
#' @aliases loadingplot,lmDME-method
stopifnot(require(pls))
setMethod(f="loadingplot",signature="lmDME", definition=function(object, term.x,
  term.y, pc=1, ord.x, col, ...){
  ##Check for PCA decomposition
  if(length(componentsType(object))==0){
     stop("No available PCA decomposed model. Please run decomposition(x)")
  }
  if(componentsType(object)["decomposition"]!="pca"){
     stop(paste("No available PCA decomposed model. Please run",
      "decomposition(x,decomposition=\"pca\")"))
  }

  ##Check term presence and appropriate order
  ## Term parameter were passed to the function
  stopifnot(!missing(term.x) & !missing(term.y))
  term <- c(paste(term.x,":",term.y,sep=""), paste(term.y,":",term.x,sep=""))
  if(!any(term %in% names(object@components))){
    stop("Term specification ", term[1], " not in modelDecomposition(object).",
      "Maybe misspelled?")}
  
  ##Get term positions and levels
  termPos<-which(names(object@components) %in% term)
  termxPos<-which(term %in% names(object@components))
  levelsY<-levels(object@design[, term.y])
  levelsX<-levels(object@design[, term.x])

  ##Get the loadings and parse them for plotting
  pcaLoad<-loadings(object@components[[termPos]])
  pcaLoad<-sapply(levelsY, function(y){
    sapply(levelsX,function(x){
      pos<-switch(termxPos, 
        '1'=rownames(pcaLoad)==paste(term.x, x, ":", term.y, y, sep=""), 
        '2'=rownames(pcaLoad)==paste(term.y, y, ":", term.x, x, sep=""))
      pcaLoad[pos,pc]
    })
  })
  
  ##Finally the loadingplot
  if(missing(ord.x)){ord.x<-1:nrow(pcaLoad)}
  if(missing(col)){col<-1:ncol(pcaLoad)}
    variance<-round(object@components[[termPos]]$sdev[pc] /
      sum(object@components[[termPos]]$sdev)*100, digits=2)
    matplot(pcaLoad[ord.x,], type="b", main=term[termxPos], axes=FALSE, 
      ylab=paste("PC-", pc, "- Exp. Var. (", variance, "%)", sep=""),
      xlab=term.x, pch=16, col=col, ...)
  box()
  axis(side=2)
  axis(1, at=1:nrow(pcaLoad), labels=rownames(pcaLoad)[ord.x])
  legend("topright", legend=paste(term.y, levelsY, sep="="), pch=16,
    box.col="transparent", col=col)
})
