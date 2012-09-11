#' \code{decomposition} of lmDME object
#' 
#' This function calculates the decomposition of variance or covariance structure using Principal Components Analysis (PCA) or Partial Least Squared Regression (PLSR), over ANOVA decomposed lmDME object. In this context, in a two factor experimental design with interaction, the linear model of the i-th observation (gene) can be written: \cr \eqn{X=\mu+XA+XB+XAB+\epsilon} \cr where \itemize{ \item X stands for the observed value \item the intercept the intercept \eqn{\mu} \item XA, XB and XAB are the first, second and interaction coefficients respectively \item The error term \eqn{\epsilon ~ N(0,\sigma^2)}.} The the model is iterative decomposed in a step by step fashion decomposing one term at each time by calling lmdme constructor: \cr \enumerate{ \item Step 1: \eqn{X=\mu+E1} \item Step 2: E1=XA+E2 \item Step 3: E2=XB+E3 \item Step 4: E3=XAB+E4.} Then, if we apply PCA over the i-th step using Ei-1 matrix it is known as \emph{APCA} but if applied over the coefficients Xi it is called \emph{ASCA}. The same decomposition schema can also be used with PLSR.
#'
#' @param object lmDME class object.
#' @param decomposition character to indicate the decomposition to be carry out, i.e., \emph{PCA} or \emph{PLSR}. Default value is PCA.
#' @param term character specifying the model term to perform the decomposition (e.g. "time" or "time:concentration" for interaction term). If term is not specified (i.e. missing) it performs the analysis over all the model terms. 
#' @param subset subset of indididuals (rows) to be included in the analysis. By default all the indididuals are included.
#' @param type character to indicate over which regressor matrix perform the decomposition ("coefficient" or "residual"). The intercept term is not included in the results, as it can be directly analyzed with the original M data.frame. Default value is "coefficient" a.k.a. ASCA. Note that "residual" performs PCA or PLS over the i-th residual Ei-1=Xi+Ei and not the residuals of the i-th model Ei.
#' @param scale character "row", "column" or "none" to indicate if the matrix should be scaled by the row, column or not  respectively. Default value is "none".
#' @param Omatrix the output matrix for PLSR only. If the parameter is missing, the output matrix  will be an identity matrix for the ASCA, otherwhise the design matrix corresponding to the specified term for APCA.
#' @param ... additional parameters for prcomp or plsr functions, according to decomposition call.
#'
#' @return Internal update of the "components" slot of the lmDME object, which is a list of \code{\link{prcomp}} or a list of mvr (\code{\link{plsr}}) objects using the given term parameter. If missing(term), the length of the list equal the number of decomposed models minus the Intercept term for coefficients or the length of decomposed models for residual decomposition. 
#'
#' @seealso \code{\link{prcomp}}, \code{\link{plsr}}
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#'
#' @references 
#' \enumerate{
#'  \item Smilde AK, Jansen JJ, Hoefsloot HCJ, Lamer RAN, Van der Greef J, Timmerman ME (2005) ANOVA-simultaneaus component analysis (ASCA): a new tool for analyzing designed metabolomics data, Bioinformatics 21,13,3043 DOI:/10.1093/bioinformatics/bti476
#'  \item Zwanenburg G, Hoefsloot HCJ, Westerhuis JA, Jansen JJ, Smilde AK (2011) ANOVA.Principal component analysis and ANOVA-simultaneaus component analysis: a comparison J. Chemometrics 25:561-567 DOI:10.1002/cem.1400
#' }
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
#' ##Just a copy of the same fit object and to perform analysis over those subjects/genes where at least one interaction coefficient is statistically different from zero (F-test over the coefficients).
#' asca <- fit; apca <-fit; plsr.residuals <- fit; plsr.coefficients <-fit
#' id<-F.p.values(fit,term="time:oxygen")[[1]]<0.001
#' ##ASCA and APCA decomposition for every available term.
#' decomposition(asca,decomposition="pca", subset=id, scale="row")
#' decomposition(apca,decomposition="pca", subset=id, scale="row", type="residual")
#' asca <- components(asca) ##get the coefficients prcomp objects 
#' apca <- components(apca) ##get the residuals prcomp objects 
#' ##PLSR decomposition for residuals and coefficients
#' decomposition(plsr.coefficients,decomposition="plsr", subset=id, scale="row")
#' decomposition(plsr.residuals,decomposition="plsr", subset=id, scale="row", type="residual")
#' plsr.coefficients <- components(plsr.coefficients) ##get the coefficients plsr objects 
#' plsr.residuals <- components(plsr.residuals) ##get the residuals plsr objects 
#' 
#' @exportMethod decomposition
#' @docType methods
#' @name decomposition
#' @rdname lmDME-decomposition
#' @aliases decomposition-methods
setGeneric(name="decomposition",def = function(object, decomposition=c("pca","plsr"), term=NULL, subset=1:nrow(object@residuals[[1]]), type = c("coefficient", "residual"), scale=c("none","row","column"), Omatrix, ...){standardGeneric("decomposition")})
#'
#' @name decomposition
#' @rdname lmDME-decomposition
#' @inheritParams decomposition
#' @aliases decomposition,lmDME-method
setMethod(f="decomposition",signature = "lmDME", definition = function(object, decomposition=c("pca","plsr"), term=NULL, subset=1:nrow(object@residuals[[1]]), type = c("coefficient", "residual"), scale=c("none","row","column"), Omatrix, ...){
  ##Obtain object name for future update
  nameObject <- deparse(substitute(object)) 
  
  ##Check parameters
  stopifnot(scale[1] %in% c("none","row","column"))
  stopifnot(type[1] %in% c("coefficient", "residual"))
  stopifnot(decomposition[1] %in% c("pca","plsr"))
  if(missing(Omatrix)){Omatrix <- NULL}
  #No sense to decompose the intercept
  if(!is.null(term)){
    if(term == "(Intercept)"){stop("Error: does not make sence to decompose the intercept coefficients/residuals. Use the centered original data object")}
    else{
       if(!any(term %in% names(object@modelDecomposition))){stop("ERROR: ", term, " not in model: ", object@model)}
    }
  }

  ##Over which matrix is to be the carry the scaletion and calculate it
  data <- switch(type[1], 
		coefficient = coefficients(object,term), 
		residual = {
			      ##Get the names of the term if missing
			      if(is.null(term)){term <- names(object@modelDecomposition)[-1]}			  
			      ##Get the residual of the previous model
			      datum <- as.list(term)
			      names(datum) <- term
			      lapply(datum,function(x){
				       index <- which(names(object@modelDecomposition)==x)
				       ##Get the previous residual
				       aux <- object@residuals[[index-1]]
				       colnames(aux) <- colnames(object@residuals[[index]])
				       aux
				    })
			     }
	  )
  
  ##Remove the intercept term if available
  if("(Intercept)" %in% names(data) & type[1] == "coefficient"){data <- data[-1]} 

  ##Scale if required
  data <- lapply(data, function(datum){
	    switch(scale[1],none = datum[subset,],row  = t(scale(t(datum[subset,]))),column = scale(datum[subset,]))
	  })

  ##decomposition calculation
  index <- as.list(1:length(data))
  names(index) <- names(data)
  object@components <- lapply(index, function(datum){
      switch(decomposition[1], pca= prcomp(data[[datum]],...),
	    plsr={
		  data.pls <- list(X=t(data[[datum]]))
		  
		  ##Check for Omatrix presence
		  if(is.null(Omatrix)){
		      Omatrix <- switch(type[1], coefficient = diag(1,nrow(data.pls$X)),
				 residual =  model.matrix(object@modelDecomposition[[names(data)[datum]]],object@design))
 		      rownames(Omatrix)   <- switch(type[1], coefficient = rownames(data.pls$X),residual = NULL)    
		      rownames(data.pls$X)<- switch(type[1], coefficient = rownames(data.pls$X),residual = NULL)
		  }
		  data.pls$Y <- Omatrix
		  plsr(formula=Y~X,data=data.pls,...)
		}#from plsr
	    )#from decomposition 
  })

  ##Update lmDME object
  object@componentsType <- c(decomposition=decomposition[1], type=type[1],scale=scale[1])
  assign(nameObject,object,envir=parent.frame())

  return(invisible(NULL))
})
