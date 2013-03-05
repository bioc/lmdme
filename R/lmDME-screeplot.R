#' Plot a \code{screeplot} of a PCA decomposed lmDME object
#' 
#' Screeplot over each decomposed "pca" model present in lmDME components slot.
#'
#' @param x lmDME class object.
#' @param independent logical indicating whether the screeplots should be
#'  plotted together. Default value is FALSE.
#' @param col which color to use for each decomposed model. Default value
#'  seq(along= components(x)).
#' @param npcs integer with the number of components to plot. By default all
#'  present components are plotted.
#' @param term character with the corresponding term/s for biploting. Default
#'  value is NULL in order to obtain every available biplot/s.
#' @param mfcol numeric vector for par layout. If missing mfcol=c(1,2) will be
#'  used if more than one biplot is available. Use  mfcol==NULL to override par
#'  call inside biplot function.
#' @param ... additional parameters for screeplot or plot/lines according to
#'  independent FALSE or TRUE respectively.
#'
#' @return plotted screeplot/s of the components slot if PCA decomposition was
#'  applied.
#'
#' @seealso stats::screeplot
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
#' design$time <-as.factor(design$time)
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
#' id<-F.p.values(fit,term="time:oxygen")[[1]]<0.001
#' decomposition(fit, decomposition="pca", scale="row", subset=id) 
#'
#' \dontrun{
#' par(mfrow=c(2,2))
#'
#' ##Do not call par inside
#' screeplot(fit,mfcol=NULL)
#'
#' ##Just the term of interest
#' screeplot(fit,term="time")
#'
#' ##In separate graphics
#' screeplot(fit,term=c("time","oxygen"),mfcol=c(1,1))
#'
#' ##All term in the same graphic device
#' screeplot(fit,mfcol=c(1,3))
#'
#' ##All in the same graphic
#' screeplot(fit,independent=FALSE)
#' }
#' }
#' 
#' @exportMethod screeplot
#' @docType methods
#' @name screeplot
#' @rdname lmDME-screeplot
#' @usage \S4method{screeplot}{lmDME}(x, independent=TRUE, col=seq(along=components(x)), npcs, term=NULL, mfcol, ...)
#' @aliases screeplot,lmDME-method
setMethod(f="screeplot", signature="lmDME", definition=function(x,
  independent=TRUE, col=seq(along=components(x)), npcs, term=NULL, mfcol, ...){
  ##Check is componentType is available
  if(length(componentsType(x))==0){
    stop("No available PCA decomposed model. Please run decomposition(x)")
  }
  if(componentsType(x)["decomposition"]!="pca"){
    stop(paste("No available PCA decomposed model. Please run",
      "decomposition(x,decomposition=\"pca\")"))
  }

  ##Get the components to use
  component<-components(x, term=term)

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

  ##Should the screeplot be plotted independent
  if(independent){
    ##Check npcs
    if(missing(npcs)){npcs<-NULL}

    invisible(sapply(seq(along=component),function(pca){  
      ##Tunning the graphic device if required
      if(!is.null(mfcol) & ((pca-1) %% (graphicsPerPar))==0){
        X11()
        par(mfcol=mfcol)
      }

      ##Check npcs
      if(is.null(npcs)){
        screeplot(component[[pca]], main=names(component)[[pca]], ...,
          col=col[pca])
      }else{
        screeplot(component[[pca]], main=names(component)[[pca]], npcs=npcs,
          col=col[pca], ...)
      }
   }))#sapply(seq(along=component)
  }else{
      ##Not independent screeplots (same device)
      variance<-lapply(component, function(pca){pca$sdev^2/sum(pca$sdev^2)})
      maxVariance<-max(unlist(lapply(variance,max)))
      maxComp<-ifelse(!missing(npcs), npcs, max(unlist(lapply(variance,
        length)))) 

      ##First screeplot
      plot(x=seq(along=variance[[1]]), y=variance[[1]], type="b",
        ylab="Variance", xlab="Component", xlim=c(1, maxComp), ylim=c(0,
        maxVariance), col=col[1], ...)

      ##The rest if available
      invisible(sapply(seq(along=variance)[-1], function(pca, variance){
        lines(x=seq(along=variance[[pca]]), y=variance[[pca]],
        type="b", col=col[pca], ...)
      },variance))
      legend("topright", legend=names(component), lty="solid", 
        box.col="transparent", col=col)
  }
})  
