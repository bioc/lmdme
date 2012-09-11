\name{fitted.values}
\alias{F.p.values}
\alias{F.p.values,lmDME-method}
\alias{coef}
\alias{coef,lmDME-method}
\alias{coefficients}
\alias{coefficients,lmDME-method}
\alias{components}
\alias{components,lmDME-method}
\alias{componentsType}
\alias{componentsType,lmDME-method}
\alias{design}
\alias{design,lmDME-method}
\alias{fitted}
\alias{fitted,lmDME-method}
\alias{fitted.values}
\alias{fitted.values,lmDME-method}
\alias{model}
\alias{model,lmDME-method}
\alias{modelDecomposition}
\alias{modelDecomposition,lmDME-method}
\alias{p.values}
\alias{p.values,lmDME-method}
\alias{resid}
\alias{resid,lmDME-method}
\alias{residuals}
\alias{residuals,lmDME-method}
\title{Getters for lmDME object}
\usage{
  \S4method{fitted.values}{lmDME}(object,term)

  \S4method{fitted}{lmDME}(object,term)

  \S4method{coef}{lmDME}(object,term)

  \S4method{coefficients}{lmDME}(object,term)

  \S4method{resid}{lmDME}(object,term)

  \S4method{residuals}{lmDME}(object,term)

  F.p.values(object, term = NULL)

  \S4method{F.p.values}{lmDME}(object, term = NULL)

  p.values(object, term = NULL)

  \S4method{p.values}{lmDME}(object, term = NULL)

  modelDecomposition(object, term = NULL)

  \S4method{modelDecomposition}{lmDME}(object, term = NULL)

  components(object, term = NULL)

  \S4method{components}{lmDME}(object, term = NULL)

  componentsType(object, term = NULL)

  \S4method{componentsType}{lmDME}(object, term = NULL)

  model(object)

  \S4method{model}{lmDME}(object)

  design(object)

  \S4method{design}{lmDME}(object)
}
\arguments{
  \item{object}{lmDME class object.}

  \item{term}{character with the corresponding term/s to
  return. Default value is NULL in order to return every
  available term/s.}
}
\value{
  according to the call one of the following objects can be
  returned \item{design}{used experiment design
  data.frame.} \item{model}{used decompose formula.}
  \item{modelDecomposition}{list of decomposed model
  formulas.} \item{residuals or resid, coef or
  coefficients, fitted or fitted.values, p.values or
  F.p.values}{list of appropieate slot where each item is a
  matrix that will have G rows(individuals) x k
  columns(levels of the corresponding model term).}
  \item{components}{list with corresponding PCA or PLSR
  term according to decomposition function call.}
  \item{componentsType}{character name vector with the
  information of the components calculation.}
}
\description{
  Obtain lmDME slot information, acording to the given
  function call (see Values). If term parameter is not
  specified, it will return all the available terms,
  otherwise just the one specified.
}
\examples{
if (require(stemHypoxia)) {
  data(stemHypoxia)
  ##Just to make a balance dataset in the Fisher sense (2 samples per time*oxygen levels)
  design<-design[design$time \%in\% c("0.5","1","5") & design$oxygen \%in\% c("1","5","21"),]
  design$time <-as.factor(design$time)
  design$oxygen<-as.factor(design$oxygen)
  rownames(M)<-M[,1]
  M <- M[,colnames(M) \%in\% design$samplename] #Keeping appropiate samples only
  ##ANOVA decomposition
  fit <- lmdme(model=~time+oxygen+time:oxygen,data=M,design=design)
  fit.model  <- model(fit) #The model formula used
  fit.design <- design(fit)#The design data.frame used
  fit.modelDecomposition<- modelDecomposition(fit) #How the decomposition was carried out
  ##Getting the coefficients
  timeCoef <- coef(fit,term="time")  #or coef(fit) to get all term coefficients
  fit.p.values <- p.values(fit,term="time") #for the same term
  fit.f.values <- F.p.values(fit,term="time") #or the F-values
  ##Getting the residuals or fitted values
  interactionResid <- resid(fit, term="time:oxygen") #or resid(fit) to get all term residuals
  oxygenDFit <- fitted(fit, term="(Intercept)") #or fitted(fit) to get all term fitted values
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{lmdme}}, \code{\link{decomposition}},
  \code{\link{print}}, \code{\link{show}}
}

