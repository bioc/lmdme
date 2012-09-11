\name{permutation}
\alias{permutation}
\alias{permutation,formula,data.frame,data.frame-method}
\alias{permutation-methods}
\title{\code{permutation} of the specified Linear Model for Designed Multivariate Experiment (lmDME)}
\usage{
  \S4method{permutation}{formula,data.frame,data.frame}(model,
    data, design, Bayes = FALSE, verbose = FALSE,
    NPermutations=100,nCpus=1, ...)
}
\arguments{
  \item{model}{formula object to carry out the
  decomposition}

  \item{data}{data.frame with individuals (rows) and
  samples/conditions (columns)}

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

  \item{NPermutations}{number of permutations to be
  calculated. Default value is 100.}

  \item{nCpus}{number of cores to be used. Default value is
  1, i.e. sequencial calculation.}
}
\value{
  \item{list}{contains the original lmDME object plus the
  required amount of permuted versions}
}
\description{
  Produces de specified lmDME plus the requied permuted
  objects (sampling the columns), usign the same parameters
  to fit the additional models.
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
  ##Just to test if it works. On real scenario use NPermutations >= 100 if the conditions (columns) of M allows it.
  permuted <- permutation(model=~time*oxygen,data=M,design=design, NPermutations=2,nCpus=3)##verbose=FALSE
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{lmdme}}
}

