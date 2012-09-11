\name{print}
\alias{print}
\alias{print,lmDME-method}
\alias{show}
\alias{show,lmDME-method}
\alias{summary}
\alias{summary,lmDME-method}
\title{\code{Show}, \code{Print}  or \code{Summary} a lmDME object}
\usage{
  \S4method{print}{lmDME}(x,term)

  \S4method{show}{lmDME}(object)

  \S4method{summary}{lmDME}(object)
}
\arguments{
  \item{x}{lmDME class object}

  \item{object}{lmDME class object}

  \item{term}{character with the corresponding term to
  return. Default value is NULL to return every decomposed
  term (if more than one available)}
}
\value{
  according to the call \item{print, show or
  summary:}{console output text with increasing detail of
  lmDME object.} \item{show or summary}{console output text
  of the lmDME object, plus a data.frame with model
  decomposition summary data.}
}
\description{
  Generic Show/Print/Summary Method for lmDME class output
  visualization.
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
}
\dontrun{
  fit #equivalent to call show(fit)
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

