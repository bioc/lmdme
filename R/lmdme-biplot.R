#' Plot a \code{biplot} of a lmdme object
#' 
#' Plot a biplot over each decomposed "pca" or "plsr" present in lmdme
#' component object's slot.
#'
#' @param x lmdme class object.
#' @param comp a two component vector with the PC components to plot. Default 
#'  comp=1:2.
#' @param xlab character for the x-label title for PCA biplots.
#' @param ylab character for the y-label title for PCA biplots.
#' @param term character with the corresponding term/s for biploting. Default 
#'  value is NULL in order to obtain every available biplot/s.
#' @param mfcol numeric vector for par layout. If missing mfcol=c(1,2) will be
#'  used if more than one biplot is available. Use  mfcol==NULL to override par
#'  call inside biplot function.
#' @param xlabs,ylabs vector of character strings to label the first/second set
#'  of points. The default is to use dimname of "x"/"y", or "1:n" if the dimname
#'  is NULL for the respective set of points. If a single character is
#'  passed e.g. "o", the same character is used for all the points.
#' @param which character to indicate the type of biplot to use when plsr
#'  decomposition is applied. Default value is "x" (X scores and loadings), "y"
#'  for (Y scores and loadings), "scores" (X and Y scores) or "loadings" (X and
#'  Y loadings). See \code{\link{biplot.mvr}} for details.
#' @param ... additional parameters for \code{\link{biplot.prcomp}}(pca) or
#'  \code{\link{biplot.mvr}}(plsr)
#'
#' @return plotted biplot/s of the component/s of the given lmdme object. If
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
#' ##Just to make a balanced dataset in the Fisher sense (2 samples per 
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
#' ##on coefficients).
#' id<-F.p.values(fit, term="time:oxygen")<0.001
#' decomposition(fit, decomposition="pca",scale="row",subset=id) 
#' 
#' \dontrun{
#' ##Does not call par inside
#' par(mfrow=c(2,2))
#' biplot(fit, xlabs="o", mfcol=NULL) 
#' 
#' ##Just the term of interest
#' biplot(fit, xlabs="o", term="time")
#'
#' ##In separate graphics
#' biplot(fit, xlabs="o", term=c("time", "oxygen"), mfcol=c(1,1))
#' 
#' ##All terms in the same graphic
#' biplot(fit, xlabs="o", mfcol=c(1,3))
#' }
#' }
#'
#' ##Now using plsr on interaction coefficients
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
#' biplot(fit, which="loadings", mfcol=NULL, xlabs="o")
#' }
#'
#' @exportMethod biplot
#' @docType methods
#' @usage \S4method{biplot}{lmdme}(x, comp=1:2, xlab=NULL, ylab=NULL, term=NULL, mfcol, xlabs, ylabs, which, ...)
#' @name biplot
#' @rdname lmdme-biplot
#' @aliases biplot,lmdme-method
setMethod(f="biplot", signature="lmdme", definition=function(x, comp=1:2,
  xlab=NULL, ylab=NULL, term=NULL, mfcol, xlabs, ylabs, which, ...){
  ##Check that is at least one biplot available
  component<-components(x, term=term, drop=FALSE)
  stopifnot(length(component)!=0)

  ##Check for xlabs, ylabs and which presence
  if(missing(xlabs)){xlabs<-NA_character_}
  if(missing(ylabs)){ylabs<-NA_character_}
  if(missing(which)){which<-"x"}
   
  ##Check mfcol to adjust par parameter and create the first biplot
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

    ##Check xlabs presence
    xData<-switch(x@componentsType["decomposition"],
      pca=component[[index]]$x,
      plsr={switch(which,
              loadings=component[[index]]$loadings,
              scores=component[[index]]$scores,
              y=component[[index]]$Yscores,
              x=component[[index]]$scores) 
    }) 
    if(!all(is.na(xlabs))){
      if(length(xlabs)==1){
        xlabs<-rep(xlabs,nrow(xData))
      }
    }else{
      ##Check default behavior 
      if(is.null(dimnames(xData))){
        xlabs<-1:nrow(xData)
      }else{
        xlabs<-rownames(xData)
      }
    }

   ##Check ylabs presence
   yData<-switch(x@componentsType["decomposition"],
      pca=component[[index]]$rotation,
      plsr={switch(which,
              loadings=component[[index]]$Yloadings,
              scores=component[[index]]$Yscores,
              y=component[[index]]$Yloadings,
              x=component[[index]]$loadings) 
    }) 
    if(!all(is.na(ylabs))){
      if(length(ylabs)==1){
        ylabs<-rep(ylabs,nrow(yData))
      }
    }else{
      ##Check default behavior 
      if(is.null(dimnames(yData))){
        ylabs<-1:nrow(yData)
      }else{
        ylabs<-rownames(yData)
      }
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
          xlab=xlab, ylab=ylab, xlabs=xlabs, ylabs=ylabs, ...)
     },
     plsr=biplot(component[[index]], comps=comp,main=names(component)[[index]],
       which=which, xlabs=xlabs, ylabs=ylabs, ...)
   )
    
    ##Just to get the reference
    abline(h=0)
    abline(v=0)
    return(NULL)
   }, xlab, ylab))
})
