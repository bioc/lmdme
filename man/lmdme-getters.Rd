\docType{methods}
\name{fitted.values}
\alias{F.p.values}
\alias{F.p.values,lmdme-method}
\alias{F.p.values-methods}
\alias{coef}
\alias{coef,lmdme-method}
\alias{coefficients}
\alias{coefficients,lmdme-method}
\alias{components}
\alias{components,lmdme-method}
\alias{components-methods}
\alias{componentsType}
\alias{componentsType,lmdme-method}
\alias{componentsType-methods}
\alias{design}
\alias{design,lmdme-method}
\alias{design-methods}
\alias{fitted}
\alias{fitted,lmdme-method}
\alias{fitted.values}
\alias{fitted.values,lmdme-method}
\alias{model}
\alias{model,lmdme-method}
\alias{model-methods}
\alias{modelDecomposition}
\alias{modelDecomposition,lmdme-method}
\alias{modelDecomposition-methods}
\alias{p.values}
\alias{p.values,lmdme-method}
\alias{p.values-methods}
\alias{resid}
\alias{resid,lmdme-method}
\alias{residuals}
\alias{residuals,lmdme-method}
\title{Getters for lmdme object}
\usage{
  \S4method{fitted.values}{lmdme}(object, term=NULL,
  drop=TRUE)

  \S4method{fitted}{lmdme}(object, term=NULL, drop=TRUE)

  \S4method{coef}{lmdme}(object, term=NULL, drop=TRUE)

  \S4method{coefficients}{lmdme}(object, term=NULL,
  drop=TRUE)

  \S4method{resid}{lmdme}(object, term=NULL, drop=TRUE)

  \S4method{residuals}{lmdme}(object, term=NULL, drop=TRUE)

  F.p.values(object, term=NULL, drop=TRUE)

  \S4method{F.p.values}{lmdme}(object, term=NULL,
  drop=TRUE)

  p.values(object, term=NULL, drop=TRUE)

  \S4method{p.values}{lmdme}(object, term=NULL, drop=TRUE)

  modelDecomposition(object, term=NULL, drop=TRUE)

  \S4method{modelDecomposition}{lmdme}(object, term=NULL,
  drop=TRUE)

  components(object, term=NULL, drop=TRUE)

  \S4method{components}{lmdme}(object, term=NULL,
  drop=TRUE)

  componentsType(object)

  \S4method{componentsType}{lmdme}(object)

  model(object)

  \S4method{model}{lmdme}(object)

  design(object)

  \S4method{design}{lmdme}(object)
}
\arguments{
  \item{object}{lmdme class object.}

  \item{term}{character with the corresponding term/s to
  return. Default value is NULL in order to return every
  available term/s.}

  \item{drop}{should try to drop list structure if
  length==1? Default value is TRUE}
}
\value{
  according to the call one of the following objects can be
  returned \item{design}{experiment design data.frame
  used.} \item{model}{decomposed formula used.}
  \item{modelDecomposition}{list of decomposed model
  formulas.} \item{residuals, resid, coef, coefficients,
  fitted, fitted.values, p.values or F.p.values}{list of
  appropriate slot where each item is a matrix that will
  have G rows (individuals) x k columns (levels of the
  corresponding model term).} \item{components}{list with
  corresponding PCA or PLSR terms according to the
  decomposition function call.}
  \item{componentsType}{character name vector with the
  information of the component calculations.}
}
\description{
  To obtain lmdme slot information, according to the given
  function call (see Values). If a term parameter is not
  specified, it will return all the available terms.
  Otherwise, just the one specified.
}
\examples{
{
data(stemHypoxia)

##Just to make a balanced dataset in the Fisher sense (2 samples per
## time*oxygen levels)
design<-design[design$time \%in\% c(0.5, 1, 5) & design$oxygen \%in\% c(1,5,21),]
design$time<-as.factor(design$time)
design$oxygen<-as.factor(design$oxygen)
rownames(M)<-M[, 1]

##Keeping appropriate samples only
M<-M[, colnames(M) \%in\% design$samplename]

##ANOVA decomposition
fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)

##Let's inspect how the decomposition process was carried out:
##a) The model formula used
##b) The design data.frame used
##c) The decomposition itself
fit.model<-model(fit)
fit.design<-design(fit)
fit.modelDecomposition<-modelDecomposition(fit)

##Getting the specific "time" term coefficients, p-values or F-values.
## Omit "term" parameter for all available terms.
timeCoef<-coef(fit,term="time")
fit.p.values<-p.values(fit,term="time")
fit.f.values<-F.p.values(fit,term="time")

##Getting the residuals or fitted values, for the interaction "time:oxygen"
## term. Omit "term" parameter for all available terms.
interactionResid<-resid(fit, term="time:oxygen")
interactionFit<-fitted(fit, term="time:oxygen")
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{lmdme}}, \code{\link{decomposition}},
  \code{\link{print}}, \code{\link{show}}
}

