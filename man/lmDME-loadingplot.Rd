\name{loadingplot}
\alias{loadingplot}
\alias{loadingplot,lmDME-method}
\title{\code{loadingplot} of interaction PCA decomposed lmDME object}
\usage{
  \S4method{loadingplot}{lmDME}(object, term.x, term.y,
    pc=1, ord.x, col, ...)
}
\arguments{
  \item{object}{lmDME class object.}

  \item{term.x,term.y}{character indicating the model
  principal factor for the interaction term (term.x:term.y
  or term.y:term.x) for the corresponding x or y axis.}

  \item{pc}{integer indicating which principal component
  loading is to be plotted on the y-axis. Default value is
  1.}

  \item{col}{which color to use for each level present in
  term.y.}

  \item{ord.x}{numeric indicating the term.x levels order,
  for plotting purpoouses. If missing the levels order is
  used.}

  \item{...}{additional parameters for maplot.}
}
\value{
  loading plot of the selected interaction (term.x:term.y)
  components lmDME object's slot, if PCA decomposition was
  applied.
}
\description{
  This function plots the PCA loadings for a given
  interaction (A:B) lmDME object's components slot, for the
  given "pc" component. The user can choose which term (A
  or B) is used for x-axis and y-axis functions (B or A)
  respectively.
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
  decomposition(fit,decomposition="pca",scale="row", subset=id)
}
\dontrun{
 loadingplot(fit, term.x="time", term.y="oxygen")
 loadingplot(fit, term.x="oxygen", term.y="time")	##Or change the axis order
 loadingplot(fit, term.x="time", term.y="oxygen",pc=2)	##Or change the PC to display
 loadingplot(fit, term.x="time", term.y="oxygen",ord.x=3:1)	##Or the order of x-levels
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}

