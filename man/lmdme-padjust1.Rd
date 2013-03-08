\docType{methods}
\name{p.adjust}
\alias{lmdme-padjust}
\alias{p.adjust}
\alias{p.adjust,ANY-method}
\alias{p.adjust,lmdme-method}
\alias{p.adjust-methods}
\title{\code{p.adjust} of p-values for Multiple Test Comparisons Correction}
\usage{
  \S4method{p.adjust}{lmdme}(p, term=NULL,
  method=p.adjust.methods, drop=TRUE)
}
\arguments{
  \item{p}{numeric vector of p-values as in stats::p.adjust
  or lmdme class object.}

  \item{method}{correction method available in
  \code{\link{p.adjust.methods}}.}

  \item{term}{character with the corresponding term to
  return.}

  \item{...}{other arguments.}

  \item{drop}{should try to drop list structure if
  length==1? Default value is TRUE}
}
\value{
  according to the call one of the following objects can be
  returned \item{numeric}{vector of adjusted p-values.}
  \item{matrix}{for lmdme object If term!=NULL, the
  corresponding character is looked up within list of
  p.values returning the associated matrix of G
  rows(individuals) x k columns(levels of the corresponding
  model term) with the adjusted p-values.}
}
\description{
  Given a set of p-values, returns p-values adjusted using
  one of several methods.
}
\examples{
{
data(stemHypoxia)

##Just to make a balance dataset in the Fisher sense (2 samples per
## time*oxygen levels)
design<-design[design$time \%in\% c(0.5, 1, 5) & design$oxygen \%in\% c(1,5,21),]
design$time<-as.factor(design$time)
design$oxygen<-as.factor(design$oxygen)
rownames(M)<-M[, 1]

##Keeping appropriate samples only
M<-M[, colnames(M) \%in\% design$samplename]

##ANOVA decomposition
fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)

##Adjust p-values only over interaction p.values using false discovery rate
## method
pInteraction<-p.values(fit, term="time:oxygen")
FDRInteraction<-p.adjust(fit, term="time:oxygen", method="fdr")
corrected<-sum(pInteraction < 0.05) - sum(FDRInteraction < 0.05)
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{p.adjust}}, \code{\link{p.adjust.methods}}
}

