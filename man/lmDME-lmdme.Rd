\name{lmdme}
\alias{lmdme}
\alias{lmdme,formula,ANY,data.frame-method}
\alias{lmdme-methods}
\title{High level constructor of lmDME class object}
\usage{
  \S4method{lmdme}{formula,ANY,data.frame}(model, data,
    design, Bayes = FALSE, verbose = FALSE, ...)
}
\arguments{
  \item{model}{formula object to carry out the
  decomposition}

  \item{data}{matrix or data.frame with individuals/genes
  (rows) and samples/conditions (columns)}

  \item{design}{data.frame with the design of the
  experiment, (rows) samples/conditions as in data columns
  and as many columns to indicate the factors presents in
  each sample.}

  \item{Bayes}{Should limma estimate empirical Bayes
  statistics, i. e., moderated t-statistics? Default value
  is FALSE.}

  \item{verbose}{Should the process progress be printed?
  Default value is FALSE.}

  \item{...}{Additional parameters for lmFit function}
}
\value{
  \item{lmDME}{lmDME class object with the corresponding
  completed slots according to the given model}
}
\description{
  Linear model ANOVA decomposition of Designed Multivariate
  Experiments based on limma lmFit implementation. For
  example in a two factor experimental design with
  interaction, the linear model of the i-th observation
  (gene) can be written: \cr \eqn{X=\mu+A+B+AB+\epsilon}
  \cr where \itemize{ \item X stands for the observed value
  \item the intercept the intercept \eqn{\mu} \item A, B
  and AB are the first, second and interaction terms
  respectively \item The error term \eqn{\epsilon ~
  N(0,\sigma^2)}.} The the model is iterative decomposed in
  a step by step fashion decomposing one term at each time:
  \cr \enumerate{ \item The intercept is estimated using
  \eqn{X=\mu+E1} \item The first factor (A) using E1=A+E2
  \item The second factor (B) using E2=B+E3 \item The
  interaction (AB) using E3=AB+E4.} For each decomposed
  step the model, residuals, coefficients, p-values and
  F-values are stored in a list container, so their
  corresponding length is equal to the number of model
  terms + 1 (the intercept).
}
\section{Note}{
  use \bold{lmdme} high level constructor for the creation
  of the class instead of directely calling its constructor
  by means of new.
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
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\references{
  \enumerate{ \item Smyth, G. K. (2005). Limma: linear
  models for microarray data. In: Bioinformatics and
  Computational Biology Solutions using R and Bioconductor.
  R. Gentleman, V. Carey, S. Dudoit, R. Irizarry, W.  Huber
  (eds), Springer, New York, pages 397--420. }
}
\seealso{
  \code{\link{decomposition}}, \code{\link{lmFit}}
}

