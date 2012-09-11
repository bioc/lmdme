\name{leverage}
\alias{leverage}
\alias{leverage,lmDME-method}
\alias{leverage-methods}
\title{\code{leverage} test of Linear model for Designed Multivariate Experiments objects}
\usage{
  \S4method{leverage}{lmDME}(object, comps = 1:2, term,
    level = 0.95)
}
\arguments{
  \item{object}{lmDME class object.}

  \item{comps}{a numeric vector indicating the PCA
  component indexes to keep. Default the first two
  components (1:2)}

  \item{term}{a character specifying the model term.}

  \item{level}{the quantile level. Default value 0.95}
}
\value{
  data.frame with the following fields
  \item{leverage}{numeric for the corresponding row
  leverage} \item{over}{logical indicating if the leverage
  > quantile(leverage,level) for the given decomposed term}
}
\description{
  This function calculates the leverage test for each
  individual using the Principal Component Analysis (comps
  function) over the coefficients of the given decomposed
  model term.
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
  ##Leverages for the first two Principal Components and q95
  leverages <- leverage(fit,term="time:oxygen")
  ##Leverages for the first three Principal Components and q99
  leverages <- leverage(fit,comps=1:3,term="time:oxygen",level=0.99)
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\references{
  Tarazona S, Prado-Lopez S, Dopazo J, Ferrer A, Conesa A,
  Variable Selection for Multifactorial Genomic Data,
  Chemometrics and Intelligent Laboratory Systems,
  110:113-122 (2012)
}
\seealso{
  \code{\link{prcomp}}, \code{\link{quantile}}
}

