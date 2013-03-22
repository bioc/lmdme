\docType{methods}
\name{print}
\alias{print}
\alias{print,lmdme-method}
\alias{show}
\alias{show,lmdme-method}
\alias{summary}
\alias{summary,lmdme-method}
\title{\code{Show}, \code{Print} or \code{Summary} a lmdme object}
\usage{
  \S4method{print}{lmdme}(x, term)

  \S4method{show}{lmdme}(object)

  \S4method{summary}{lmdme}(object)
}
\arguments{
  \item{x}{lmdme class object.}

  \item{object}{lmdme class object.}

  \item{term}{character with the corresponding term to
  return. Default value is NULL to return every decomposed
  term (if more than one is available).}
}
\value{
  according to the call \item{print, show or
  summary:}{console output text with increasing detail of
  lmdme object.} \item{show or summary}{console output text
  of the lmdme object, plus a data.frame with model
  decomposition summary data.}
}
\description{
  Generic Show/Print/Summary method for lmdme class output
  visualization.
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
}
\dontrun{
#equivalent to call show(fit)
fit
print(fit)
summary(fit)
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{lmdme}}, \code{\link{coef}},
  \code{\link{resid}}, \code{\link{fitted}},
  \code{\link{modelDecomposition}},
  \code{\link{components}}, \code{\link{componentsType}}
}

