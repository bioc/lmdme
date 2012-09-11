#' High level constructor of lmDME class object
#'
#' Linear model ANOVA decomposition of Designed Multivariate Experiments based on limma lmFit implementation. For example in a two factor experimental design with interaction, the linear model of the i-th observation (gene) can be written: \cr \eqn{X=\mu+A+B+AB+\epsilon} \cr where \itemize{ \item X stands for the observed value \item the intercept the intercept \eqn{\mu} \item A, B and AB are the first, second and interaction terms respectively \item The error term \eqn{\epsilon ~ N(0,\sigma^2)}.} The the model is iterative decomposed in a step by step fashion decomposing one term at each time: \cr \enumerate{ \item The intercept is estimated using \eqn{X=\mu+E1} \item The first factor (A) using E1=A+E2 \item The second factor (B) using E2=B+E3 \item The interaction (AB) using E3=AB+E4.} For each decomposed step the model, residuals, coefficients, p-values and F-values are stored in a list container, so their corresponding length is equal to the number of model terms + 1 (the intercept). 
#' 
#' @param model formula object to carry out the decomposition
#' @param data matrix or data.frame with individuals/genes (rows) and samples/conditions (columns)
#' @param design data.frame with the design of the experiment, (rows) samples/conditions as in data columns and as many columns to indicate the factors presents in each sample.
#' @param Bayes Should limma estimate empirical Bayes statistics, i. e., moderated t-statistics? Default value is FALSE.
#' @param verbose Should the process progress be printed? Default value is FALSE.
#' @param ... Additional parameters for lmFit function
#'
#' @return 
#'  \item{lmDME}{lmDME class object with the corresponding completed slots according to the given model} 
#'   
#' @section Note: use \bold{lmdme} high level constructor for the creation of the class instead of directely calling its constructor by means of new.
#'
#' @seealso \code{\link{decomposition}}, \code{\link{lmFit}}
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#'
#' @references 
#' \enumerate{
#'  \item Smyth, G. K. (2005). Limma: linear models for microarray data. In: Bioinformatics and Computational Biology Solutions using R and Bioconductor. R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.  Huber (eds), Springer, New York, pages 397--420.
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
#'
#' @exportMethod lmdme
#' @docType methods
#' @name lmdme
#' @rdname lmDME-lmdme
#' @aliases lmdme-methods
setGeneric(name="lmdme",def = function(model,data,design,Bayes=FALSE,verbose=FALSE,...){standardGeneric("lmdme")})
#'
#' @name lmdme
#' @rdname lmDME-lmdme
#' @inheritParams lmdme
#' @aliases lmdme,formula,ANY,data.frame-method
setMethod(f="lmdme",signature=signature(model="formula",data="ANY",design="data.frame"), 
  definition = function(model,data,design,Bayes=FALSE,verbose=FALSE,...){
  ## Auxiliary functions 
  ## Print if verbose==TRUE
  if(verbose){printnow <-function(...){cat(...);flush.console()}}else{printnow<-function(...){invisible(NULL)}}

  ##Obtain the model parts
  variables <- attr(terms(model),"variables")
  response  <- rownames(attr(terms(model),"factors"))[attr(terms(model),"response")]
  terminos  <- paste(response,"~",c(1,paste(attr(terms(model),"term.labels"),"-1")))

  ##Initialization of different slots
  .Object <- new("lmDME")
  .Object@design <- design
  .Object@model  <- model
  .Object@residuals   <- vector("list",length(terminos)); 
  names(.Object@residuals) <-  c("(Intercept)",attr(terms(model),"term.labels"))
  .Object@p.values    <- .Object@residuals
  .Object@F.p.values  <- .Object@residuals 
  .Object@coefficients<- .Object@residuals
  .Object@modelDecomposition <- lapply(terminos,as.formula)
  names(.Object@modelDecomposition) <- names(.Object@residuals)
  validObject(.Object)
  
  ##Removing the corresponding term in each step
  invisible(sapply(1:length(terminos),function(term){
    printnow("testing: ",terminos[term],"\n")

    ##Fitting the model
    mm <- model.matrix(object=as.formula(terminos[term]),data=design)
    fit<-lmFit(object=as.matrix(data), design = mm, ...)
    if(Bayes){
      fit <- eBayes(fit)
    }else{
    ##Without Bayes
      fit$t<-fit$coefficients/(fit$stdev.unscaled*fit$sigma)
      fit$df.prior <- 0                        

      fit$p.value <- 2*pt(-abs(fit$t),fit$df.residual)      
      F.stat <- classifyTestsF(fit, fstat.only = TRUE)
      
      fit$F <- as.vector(F.stat)
      df1 <- attr(F.stat, "df1")
      df2 <- attr(F.stat, "df2")
      if (df2[1] > 1e+06){ fit$F.p.value <- pchisq(df1 * fit$F, df1, lower.tail = FALSE)}
      else{fit$F.p.value <- pf(fit$F, df1, df2, lower.tail = FALSE)}
    } 

    ##Updating object
    .Object@residuals[[term]]   <<- data-t(apply(as.matrix(fit$coefficients),MARGIN=1,FUN=function(x,mm){mm%*%x},mm))
    .Object@coefficients[[term]]<<- fit$coefficients
    .Object@p.values[[term]]    <<- fit$p.value
    .Object@F.p.values[[term]]  <<- fit$F.p.value
    data <<- .Object@residuals[[term]]

    ##Setting names to the residuals
    if(term!=1){
      ##Build Interaction Factor to name the residuals
      form <- as.formula(terminos[term])
      factors <- rownames(attr(terms(form),"factors"))
      if(attr(terms(form),"response")>0){factors<-factors[-1]}
      ifactor<-as.factor(paste(factors[1],as.character(design[,factors[1]]),sep=""))
      i<-2
      while(i <= length(factors)){
	ifactor<- ifactor:as.factor(paste(factors[i],as.character(design[,factors[i]]),sep=""))
	i<-i+1
      }
      colnames(.Object@residuals[[term]])<<-as.character(ifactor)
    }
    return(NULL)
  }))
  
  validObject(.Object)
  return(.Object)
})
