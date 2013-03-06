#' Plot a \code{biplot} of a lmDME object
#' 
#' Plot a biplot over each decomposed "pca" or "plsr" present in lmDME
#' components object's slot.
#'
#' @param x lmDME class object.
#' @param comp a two component vector with the PC components to plot. Default 
#'  comp=1:2.
#' @param xlab character for the x-label title for PCA biplots.
#' @param ylab character for the y-label title for PCA biplots.
#' @param term character with the corresponding term/s for biploting. Default 
#'  value is NULL in order to obtain every available biplot/s.
#' @param mfcol numeric vector for par layout. If missing mfcol=c(1,2) will be
#'  used if more than one biplot is available. Use  mfcol==NULL to override par
#'  call inside biplot function.
#' @param ... additional parameters for \code{\link{biplot.prcomp}}(pca) or
#'  \code{\link{biplot.mvr}}(plsr)
#'
#' @return plotted biplot/s of the components slot of the given lmDME object. If
#'  \code{\link{par}}() is called before this function, the biplots can be
#'  arranged in the same window  
#'
#' @seealso \code{\link{prcomp}}, \code{\link{plsr}},
#'  \code{\link{biplot.princomp}}, \code{\link{biplot.mvr}}
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
#' id<-F.p.values(fit, term="time:oxygen")<0.001
#' decomposition(fit, decomposition="pca",scale="row",subset=id) 
#' 
#' \dontrun{
#' ##Do not call par inside
#' par(mfrow=c(2,2))
#' biplot(fit, xlabs=rep("o", sum(id)), mfcol=NULL) 
#' 
#' ##Just the term of interest
#' biplot(fit, xlabs=rep("o", sum(id)), term="time")
#'
#' ##In separate graphics
#' biplot(fit, xlabs=rep("o", sum(id)), term=c("time", "oxygen"), mfcol=c(1,1))
#' 
#' ##All term in the same graphic
#' biplot(fit, xlabs=rep("o", sum(id)), mfcol=c(1,3))
#' }
#' }
#'
#' ##Now using plsr over interaction coefficients
#' decomposition(fit, decomposition="plsr", term="time:oxygen", scale="row",
#'  subset=id)
#'
#' \dontrun{
#' par(mfrow=c(2,2))
#' 
#' ##plsr biplot by default which="x"
#' biplot(fit, which="x", mfcol=NULL) 
#' 
#' ##Other alternatives to which
#' biplot(fit, which="y", mfcol=NULL)
#' biplot(fit, which="scores", mfcol=NULL)
#' biplot(fit, which="loadings", mfcol=NULL, xlabs=rep("o", sum(id)))
#' }
#'
#' @exportMethod biplot
#' @docType methods
#' @usage \S4method{biplot}{lmDME}(x, comp=1:2, xlab=NULL, ylab=NULL, term=NULL, mfcol, ...)
#' @name biplot
#' @rdname lmDME-biplot
#' @aliases biplot,lmDME-method
setMethod(f="biplot", signature="lmDME", definition=function(x, comp=1:2,
  xlab=NULL, ylab=NULL, term=NULL, mfcol, ...){
  ##Check that is at least one biplot available
  component<-components(x, term=term, drop=FALSE)
  stopifnot(length(component)!=0)

  ##Check mfcol to adjust par parameter and creat the first biplot
  if(missing(mfcol)){
    if(length(component)==1){
      mfcol<-c(1,1)
    } 
    else{ 
      mfcol<-c(1,2)
    }
  }  
  graphicsPerPar<-prod(mfcol)
  
  ##One plot for each decomposed component
  invisible(sapply(1:length(component), function(index, xlab, ylab){
    ##Tunning the graphic device if required
    if(!is.null(mfcol) & ((index-1) %% (graphicsPerPar)) == 0){
      X11()
      par(mfcol=mfcol)
    }

    ##Selecting the type of biplot
    switch(x@componentsType["decomposition"],
      pca={
        variance<-round(100*component[[index]]$sdev /
        sum(component[[index]]$sdev), digits=2)
      
        ##Adjusting labels
        xlab <-ifelse(is.null(xlab), 
          paste("PC", comp[1], "(", variance[comp[1]], "%)", sep=""), xlab)
        ylab <- ifelse(is.null(ylab),
          paste("PC", comp[2], "(", variance[comp[2]], "%)", sep=""), ylab)

        biplot(component[[index]], main=names(component)[[index]], choices=comp,
          xlab=xlab, ylab=ylab, ...)
     },
     plsr=biplot(component[[index]], comps=comp, main=names(component)[[index]],
      ...)
   )
    
    ##Just to get the reference
    abline(h=0)
    abline(v=0)
    return(NULL)
   }, xlab, ylab))
})
