\name{lmDME-class}
\alias{lmDME-class}
\title{lmDME S4 class: Linear model decomposition for Designed Multivariate Experiments.}
\description{
  Linear Model ANOVA decomposition of Designed Multivariate
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
\section{Features}{
  \enumerate{ \item Flexible formula type interface, \item
  Fast limma based implementation based on lmFit, \item p
  and F-values over deflacted coefficients, \item ploting
  functions for PCA and PLS. }
}

\section{Slots}{
  \itemize{ \item design: data.frame with experimental
  design.  \item model: formula with the designed model to
  be decomposed.  \item modelDecomposition: list with the
  model formula obtained for each deflation step.  \item
  residuals: list of residual matrices G rows(genes) x N
  columns(arrays-designed meassurements).  \item
  coefficients: list of coefficient matrices. Each matrix
  will have G rows(genes) x k columns(levels of the
  corresponding model term).  \item p.values: list of
  p-value matrices.  \item F.p.values: list with
  corresponding F-p-values vectors for each individual.
  \item components: list with corresponding PCA or PLS
  components for the selected term/s.  \item
  componentsType: name character vector to keep process
  trace of the variance/covariance components slot:
  decomposition ("pca" or "pls"), type ("apca" for
  ANOVA-PCA or "asca" for ANOVA-SCA) and scale ("none",
  "row" or "column") }
}

\section{lmDME-general-functions}{
  \describe{ \item{print, show}{Basic output for lmDME
  class} \item{summary}{Basic statistics for lmDME class}
  \item{design, model, modelDecomposition, residuals,
  coefficients, p.values, F.p.values, components and
  componentsType}{Getters for their respective slots.} }
}

\section{ANOVA-linear-decomposition-functions}{
  \describe{ \item{lmdme}{Function that produce the
  complete ANOVA decomposition based on model specification
  through formula interface. Technicaly is a high level
  wrapper of initialize function.}
  \item{modelDecomposition}{Getter for decomposed used
  formula in each step} \item{p.adjust}{Adjust coefficients
  p-values for Multiple Comparisons Test.} \item{Fpvalues,
  pvalues}{Getters for corresponding model decomposed
  asociated coefficient statistics in each step, for each
  observation} \item{residuals, resid, coef, coefficients,
  fitted.values, fitted}{Getters for corresponding model
  decomposed in each step} \item{permutation}{Produces de
  specified lmDME in addition to the requied permuted
  objects (sampling the columns of data), usign the same
  parameters to fit the model.} }
}

\section{variance-covariance-decomposition-functions}{
  \describe{ \item{decomposition}{Function to perform PCA
  or PLS over the ANOVA deflacted terms. PCA can be
  performed over E1, E2 or E3 and it is refered as
  ANOVA-PCA (APCA) but, if it is performed over the
  coefficients it is refered as ANOVA-SCA (ASCA). On the
  other hand PLSR is based on pls library and if it is
  performed on coefficients (ASCA like) it uses the
  identity matrix for output co-variance maximization or
  can be carried out over the E1, E2 or E3 (APCA like)
  using the design matrix as output.}
  \item{components}{Getter for PCA or PLS decomposed
  models.} \item{componentsType}{Getter for componentsType
  slot.} \item{leverage}{Leverage calculation over PCA
  (APCA or ASCA) terms.} \item{biplot}{Biplots for PCA or
  PLSR decomposed terms.} \item{screeplot}{Screeplot over
  each decomposed PCA decomposed term.}
  \item{loadingplot}{Loafingplot for PCA interaction
  terms.} }
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\references{
  \enumerate{ \item Smilde AK, Jansen JJ, Hoefsloot HCJ,
  Lamer RAN, Van der Greef J, Timmerman ME (2005)
  ANOVA-simultaneaus component analysis (ASCA): a new tool
  for analyzing designed metabolomics data, Bioinformatics
  21,13,3043 DOI:/10.1093/bioinformatics/bti476 \item
  Zwanenburg G, Hoefsloot HCJ, Westerhuis JA, Jansen JJ,
  Smilde AK (2011) ANOVA.Principal component analysis and
  ANOVA-simultaneaus component analysis: a comparison J.
  Chemometrivcs 25:561-567 DOI:10.1002/cem.1400 \item
  Tarazona S, Prado-Lopez S, Dopazo J, Ferrer A, Conesa A
  (2012) Variable Selection for Multifactorial Genomic
  Data, Chemometrics and Intelligent Laboratory Systems,
  110:113-122 }
}
\seealso{
  \code{\link{lmdme}}, \code{\link{decomposition}},
  \code{\link{biplot}}, \code{\link{loadingplot}} and
  additional related lmDME class functions.
}

