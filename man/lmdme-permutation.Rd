\docType{methods}
\name{permutation}
\alias{permutation}
\alias{permutation,formula,data.frame,data.frame-method}
\alias{permutation-methods}
\title{\code{permutation} of the specified lmdme object}
\usage{
  \S4method{permutation}{formula,data.frame,data.frame}(model,data,design,Bayes=FALSE,verbose=FALSE,NPermutations=100,nCpus=1,...)
}
\arguments{
  \item{model}{formula object to carry out the
  decomposition.}

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

  \item{NPermutations}{number of permutations to be
  calculated. Default value is 100.}

  \item{nCpus}{number of cores to be used. Default value is
  1, i.e. sequential calculation.}

  \item{...}{Additional parameters for \code{\link{lmFit}}
  function.}
}
\value{
  \item{list}{contains the original lmdme object plus the
  required amount of permuted versions.}
}
\description{
  Produces de specified lmdme plus the required permuted
  objects (sampling the columns), using the same parameters
  to fit the additional models.
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

##Just to test if it works. On real scenario use NPermutations >= 100 if the
## conditions (columns) of M allows it. Verbose parameter is FALSE by default
permuted<-permutation(model=~time*oxygen, data=M, design=design,
NPermutations=2, nCpus=3)
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\seealso{
  \code{\link{lmdme}}
}

