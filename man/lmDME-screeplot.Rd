\docType{methods}
\name{screeplot}
\alias{screeplot}
\alias{screeplot,lmDME-method}
\title{Plot a \code{screeplot} of a PCA decomposed lmDME object}
\usage{
  \S4method{screeplot}{lmDME}(x, independent=TRUE,
  col=seq(along=components(x)), npcs, term=NULL, mfcol,
  ...)
}
\arguments{
  \item{x}{lmDME class object.}

  \item{independent}{logical indicating whether the
  screeplots should be plotted together. Default value is
  FALSE.}

  \item{col}{which color to use for each decomposed model.
  Default value seq(along= components(x)).}

  \item{npcs}{integer with the number of components to
  plot. By default all present components are plotted.}

  \item{term}{character with the corresponding term/s for
  biploting. Default value is NULL in order to obtain every
  available biplot/s.}

  \item{mfcol}{numeric vector for par layout. If missing
  mfcol=c(1,2) will be used if more than one biplot is
  available. Use mfcol==NULL to override par call inside
  biplot function.}

  \item{...}{additional parameters for screeplot or
  plot/lines according to independent FALSE or TRUE
  respectively.}
}
\value{
  plotted screeplot/s of the components slot if PCA
  decomposition was applied.
}
\description{
  Screeplot over each decomposed "pca" model present in
  lmDME components slot.
}
\examples{
{
data(stemHypoxia)

##Just to make a balance dataset in the Fisher sense (2 samples per
## time*oxygen levels)
design<-design[design$time \%in\% c(0.5,1,5) & design$oxygen \%in\% c(1,5,21), ]
design$time <-as.factor(design$time)
design$oxygen<-as.factor(design$oxygen)
rownames(M)<-M[, 1]

#Keeping appropriate samples only
M<-M[, colnames(M) \%in\% design$samplename]

##ANOVA decomposition
fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)

##ASCA for all the available terms, over those subjects/genes where at least
##one interaction coefficient is statistically different from zero (F-test
##over the coefficients).
id<-F.p.values(fit,term="time:oxygen")<0.001
decomposition(fit, decomposition="pca", scale="row", subset=id)

\dontrun{
par(mfrow=c(2,2))

##Do not call par inside
screeplot(fit,mfcol=NULL)

##Just the term of interest
screeplot(fit,term="time")

##In separate graphics
screeplot(fit,term=c("time","oxygen"),mfcol=c(1,1))

##All term in the same graphic device
screeplot(fit,mfcol=c(1,3))

##All in the same graphic
screeplot(fit,independent=FALSE)
}
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  stats::screeplot
}

