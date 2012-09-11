\name{decomposition}
\alias{decomposition}
\alias{decomposition,lmDME-method}
\alias{decomposition-methods}
\title{\code{decomposition} of lmDME object}
\usage{
  \S4method{decomposition}{lmDME}(object,
    decomposition = c("pca", "plsr"), term = NULL,
    subset = 1:nrow(object@residuals[[1]]),
    type = c("coefficient", "residual"),
    scale = c("none", "row", "column"), Omatrix, ...)
}
\arguments{
  \item{object}{lmDME class object.}

  \item{decomposition}{character to indicate the
  decomposition to be carry out, i.e., \emph{PCA} or
  \emph{PLSR}. Default value is PCA.}

  \item{term}{character specifying the model term to
  perform the decomposition (e.g. "time" or
  "time:concentration" for interaction term). If term is
  not specified (i.e. missing) it performs the analysis
  over all the model terms.}

  \item{subset}{subset of indididuals (rows) to be included
  in the analysis. By default all the indididuals are
  included.}

  \item{type}{character to indicate over which regressor
  matrix perform the decomposition ("coefficient" or
  "residual"). The intercept term is not included in the
  results, as it can be directly analyzed with the original
  M data.frame. Default value is "coefficient" a.k.a. ASCA.
  Note that "residual" performs PCA or PLS over the i-th
  residual Ei-1=Xi+Ei and not the residuals of the i-th
  model Ei.}

  \item{scale}{character "row", "column" or "none" to
  indicate if the matrix should be scaled by the row,
  column or not respectively. Default value is "none".}

  \item{Omatrix}{the output matrix for PLSR only. If the
  parameter is missing, the output matrix will be an
  identity matrix for the ASCA, otherwhise the design
  matrix corresponding to the specified term for APCA.}

  \item{...}{additional parameters for prcomp or plsr
  functions, according to decomposition call.}
}
\value{
  Internal update of the "components" slot of the lmDME
  object, which is a list of \code{\link{prcomp}} or a list
  of mvr (\code{\link{plsr}}) objects using the given term
  parameter. If missing(term), the length of the list equal
  the number of decomposed models minus the Intercept term
  for coefficients or the length of decomposed models for
  residual decomposition.
}
\description{
  This function calculates the decomposition of variance or
  covariance structure using Principal Components Analysis
  (PCA) or Partial Least Squared Regression (PLSR), over
  ANOVA decomposed lmDME object. In this context, in a two
  factor experimental design with interaction, the linear
  model of the i-th observation (gene) can be written: \cr
  \eqn{X=\mu+XA+XB+XAB+\epsilon} \cr where \itemize{ \item
  X stands for the observed value \item the intercept the
  intercept \eqn{\mu} \item XA, XB and XAB are the first,
  second and interaction coefficients respectively \item
  The error term \eqn{\epsilon ~ N(0,\sigma^2)}.} The the
  model is iterative decomposed in a step by step fashion
  decomposing one term at each time by calling lmdme
  constructor: \cr \enumerate{ \item Step 1: \eqn{X=\mu+E1}
  \item Step 2: E1=XA+E2 \item Step 3: E2=XB+E3 \item Step
  4: E3=XAB+E4.} Then, if we apply PCA over the i-th step
  using Ei-1 matrix it is known as \emph{APCA} but if
  applied over the coefficients Xi it is called
  \emph{ASCA}. The same decomposition schema can also be
  used with PLSR.
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
  ##Just a copy of the same fit object and to perform analysis over those subjects/genes where at least one interaction coefficient is statistically different from zero (F-test over the coefficients).
  asca <- fit; apca <-fit; plsr.residuals <- fit; plsr.coefficients <-fit
  id<-F.p.values(fit,term="time:oxygen")[[1]]<0.001
  ##ASCA and APCA decomposition for every available term.
  decomposition(asca,decomposition="pca", subset=id, scale="row")
  decomposition(apca,decomposition="pca", subset=id, scale="row", type="residual")
  asca <- components(asca) ##get the coefficients prcomp objects
  apca <- components(apca) ##get the residuals prcomp objects
  ##PLSR decomposition for residuals and coefficients
  decomposition(plsr.coefficients,decomposition="plsr", subset=id, scale="row")
  decomposition(plsr.residuals,decomposition="plsr", subset=id, scale="row", type="residual")
  plsr.coefficients <- components(plsr.coefficients) ##get the coefficients plsr objects
  plsr.residuals <- components(plsr.residuals) ##get the residuals plsr objects
}
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
  Chemometrics 25:561-567 DOI:10.1002/cem.1400 }
}
\seealso{
  \code{\link{prcomp}}, \code{\link{plsr}}
}

