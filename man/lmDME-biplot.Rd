\name{biplot}
\alias{biplot}
\alias{biplot,lmDME-method}
\title{Plot a \code{biplot} of a lmDME object}
\usage{
  \S4method{biplot}{lmDME}(x, comp=1:2, xlab = NULL, ylab =
    NULL, term=NULL, mfcol, ...)
}
\arguments{
  \item{x}{lmDME class object.}

  \item{comp}{a tow component vector with the PC components
  to plot. Default comp=1:2,}

  \item{xlab}{character for the x-label title for PCA
  biplots}

  \item{ylab}{character for the y-label title for PCA
  biplots}

  \item{term}{character with the corresponding term/s for
  biploting. Default value is NULL in order to obtain every
  available biplot/s.}

  \item{mfcol}{numeric vector for par layout. If missing
  mfcol=c(1,2) will be used if more than one biplot is
  available. Use mfcol==NULL to override par call inside
  biplot function.}

  \item{...}{additional parameters for biplot.prcomp (pca)
  or biplot.mrv (plsr)}
}
\value{
  ploted biplot/s of the components slot of the given lmDME
  object. If par() is called before this function, the
  biplots can be arranged in the same window
}
\description{
  Plot a biplot over each decomposed "pca" or "plsr"
  present in components lmDME object's slot.
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
  ##ASCA for all the available terms, over those subjects/genes where at least one interaction coefficient is statistically different from zero (F-test over the coefficients).
  id<-F.p.values(fit,term="time:oxygen")[[1]]<0.001
  decomposition(fit,decomposition="pca",scale="row",subset=id)
  \dontrun{
  par(mfrow=c(2,2)); biplot(fit,xlabs=rep("o",sum(id)),mfcol=NULL) ##Do not call par inside
  biplot(fit,xlabs=rep("o",sum(id)),term="time")##Just the term of interest
  biplot(fit,xlabs=rep("o",sum(id)),term=c("time","oxygen"),mfcol=c(1,1))##In separate graphics
  biplot(fit,xlabs=rep("o",sum(id)),mfcol=c(1,3))##All term in the same graphic
  }
  ##Now using plsr over interaction coefficients
  decomposition(fit,decomposition = "plsr", term ="time:oxygen",scale="row", subset=id)
}
\dontrun{
par(mfrow=c(2,2))
biplot(fit, which = "x", mfcol=NULL) ##plsr biplot by default which="x"
biplot(fit, which = "y", mfcol=NULL) ##Other alternatives to which
biplot(fit, which = "scores", mfcol=NULL)
biplot(fit, which = "loadings", mfcol=NULL, xlabs=rep("o",sum(id)))
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{prcomp}}, \code{\link{plsr}},
  \code{\link{biplot.princomp}}, \code{\link{biplot.mvr}}
}

