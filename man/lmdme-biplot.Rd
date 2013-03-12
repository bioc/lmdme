\docType{methods}
\name{biplot}
\alias{biplot}
\alias{biplot,lmdme-method}
\title{Plot a \code{biplot} of a lmdme object}
\usage{
  \S4method{biplot}{lmdme}(x, comp=1:2, xlab=NULL,
  ylab=NULL, term=NULL, mfcol, xlabs, ylabs, which, ...)
}
\arguments{
  \item{x}{lmdme class object.}

  \item{comp}{a two component vector with the PC components
  to plot. Default comp=1:2.}

  \item{xlab}{character for the x-label title for PCA
  biplots.}

  \item{ylab}{character for the y-label title for PCA
  biplots.}

  \item{term}{character with the corresponding term/s for
  biploting. Default value is NULL in order to obtain every
  available biplot/s.}

  \item{mfcol}{numeric vector for par layout. If missing
  mfcol=c(1,2) will be used if more than one biplot is
  available. Use mfcol==NULL to override par call inside
  biplot function.}

  \item{xlabs,ylabs}{vector of character strings to label
  the first/second set of points. The default is to use
  dimname of "x"/"y", or "1:n" if the dimname is NULL for
  the respectively set of points. If a single character is
  passed e.g. "o", the same character one is used for all
  the points.}

  \item{which}{character to indicate the type of biplot to
  use when plsr decomposition is applied. Default value is
  "x" (X scores and loadings), "y" for (Y scores and
  loadings), "scores" (X and Y scores) or "loadings" (X and
  Y loadings). See \code{\link{biplot.mvr}} for details.}

  \item{...}{additional parameters for
  \code{\link{biplot.prcomp}}(pca) or
  \code{\link{biplot.mvr}}(plsr)}
}
\value{
  plotted biplot/s of the components slot of the given
  lmdme object. If \code{\link{par}}() is called before
  this function, the biplots can be arranged in the same
  window
}
\description{
  Plot a biplot over each decomposed "pca" or "plsr"
  present in lmdme components object's slot.
}
\examples{
{
data(stemHypoxia)

##Just to make a balance dataset in the Fisher sense (2 samples per
## time*oxygen levels)
design<-design[design$time \%in\% c(0.5,1,5) & design$oxygen \%in\% c(1,5,21), ]
design$time<-as.factor(design$time)
design$oxygen<-as.factor(design$oxygen)
rownames(M)<-M[, 1]

#Keeping appropriate samples only
M<-M[, colnames(M) \%in\% design$samplename]

##ANOVA decomposition
fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)

##ASCA for all the available terms, over those subjects/genes where at least
##one interaction coefficient is statistically different from zero (F-test
##over the coefficients).
id<-F.p.values(fit, term="time:oxygen")<0.001
decomposition(fit, decomposition="pca",scale="row",subset=id)

\dontrun{
##Do not call par inside
par(mfrow=c(2,2))
biplot(fit, xlabs="o", mfcol=NULL)

##Just the term of interest
biplot(fit, xlabs="o", term="time")

##In separate graphics
biplot(fit, xlabs="o", term=c("time", "oxygen"), mfcol=c(1,1))

##All term in the same graphic
biplot(fit, xlabs="o", mfcol=c(1,3))
}
}

##Now using plsr over interaction coefficients
decomposition(fit, decomposition="plsr", term="time:oxygen", scale="row",
subset=id)

\dontrun{
par(mfrow=c(2,2))

##plsr biplot by default which="x"
biplot(fit, which="x", mfcol=NULL)

##Other alternatives to which
biplot(fit, which="y", mfcol=NULL)
biplot(fit, which="scores", mfcol=NULL)
biplot(fit, which="loadings", mfcol=NULL, xlabs="o")
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{prcomp}}, \code{\link{plsr}},
  \code{\link{biplot.princomp}}, \code{\link{biplot.mvr}}
}

