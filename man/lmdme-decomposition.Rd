\docType{methods}
\name{decomposition}
\alias{decomposition}
\alias{decomposition,lmdme-method}
\alias{decomposition-methods}
\title{\code{decomposition} of lmdme object}
\arguments{
  \item{object}{lmdme class object.}

  \item{decomposition}{character to indicate the
  decomposition to be carried out, i.e., \emph{"pca"} or
  \emph{"plsr"}. Default value is "pca".}

  \item{term}{character specifying the model term to
  perform the decomposition (e.g. "time" or
  "time:concentration" for interaction term). If the term
  is not specified (i.e. missing) it performs the analysis
  over all the model terms.}

  \item{subset}{subset of individuals (rows) to be included
  in the analysis. By default all the individuals are
  included.}

  \item{type}{character to indicate on which regression
  matrix ("coefficient" or "residual") the decomposition
  will be performed. The intercept term is not included in
  the results, as it can be directly analyzed with the
  original M data.frame. Default value is "coefficient"
  a.k.a. ASCA. Note that "residual" performs PCA or PLS on
  the i-th residual \eqn{E_{i-1}=X_{i}+E_{i}} and not the
  residuals of the i-th model \eqn{E_{i}}.}

  \item{scale}{character "row", "column" or "none" to
  indicate if the matrix should be scaled by the row,
  column or not respectively. Default value is "none".}

  \item{Omatrix}{the output matrix for PLSR only. If the
  parameter is missing, the output matrix will be an
  identity matrix for the ASCA. Otherwise, is the design
  matrix corresponding to the specified term for APCA.}

  \item{...}{additional parameters for \code{\link{prcomp}}
  or \code{\link{plsr}} functions, according to the
  decomposition call.}
}
\value{
  Internal update of the "components" slot of the lmdme
  object, which is a list of \code{\link{prcomp}} or a list
  of mvr (\code{\link{plsr}}) objects using the given term
  parameter. If missing(term), the length of the list
  equals the number of decomposed models minus the
  Intercept term for coefficients or the length of
  decomposed models for residual decomposition.
}
\description{
  This function calculates the decomposition of variance or
  covariance structure using Principal Component Analysis
  (PCA) or Partial Least Squared Regression (PLSR), on the
  ANOVA decomposed lmdme object. In this context, in a two
  factor experimental design with interaction, the linear
  model of the i-th observation (gene) can be written: \cr
  \eqn{X=\mu+X_{A}+X_{B}+X_{AB}+\epsilon} \cr where
  \itemize{ \item X stands for the observed value. \item
  the intercept \eqn{\mu}. \item \eqn{X_{A}}, \eqn{X_{B}}
  and \eqn{X_{AB}} are the first, second and interaction
  coefficients respectively. \item The error term
  \eqn{\epsilon ~ N(0,\sigma^2)}. } The model is
  iteratively decomposed in a step by step fashion,
  decomposing one term each time by calling
  \code{\link{lmdme}} constructor: \enumerate{ \item Step
  1: \eqn{X=\mu+E_{1}} \item Step 2:
  \eqn{E_{1}=X_{A}+E_{2}} \item Step 3:
  \eqn{E_{2}=X_{B}+E_{3}} \item Step 4:
  \eqn{E_{3}=X_{AB}+E_{4}} } Then, if we apply PCA on the
  i-th step using \eqn{E_{i-1}} matrix it is known as
  \emph{APCA} but if applied on the coefficients
  \eqn{X_{i}} it is called \emph{ASCA}. The same
  decomposition schema can also be used with PLSR.
}
\examples{
{
data(stemHypoxia)

##Just to make a balanced dataset in the Fisher sense (2 samples per
## time*oxygen levels)
design<-design[design$time \%in\% c(0.5,1,5) & design$oxygen \%in\% c(1,5,21), ]
design$time<-as.factor(design$time)
design$oxygen<-as.factor(design$oxygen)
rownames(M)<-M[, 1]

#Keeping appropriate samples only
M<-M[, colnames(M) \%in\% design$samplename]

##ANOVA decomposition
fit<-lmdme(model=~time+oxygen+time:oxygen, data=M, design=design)

##Just a copy of the same fit object and to perform analysis on those
##subjects/genes where at least one interaction coefficient is statistically
##different from zero (F-test on the coefficients).
asca<-fit
apca<-fit
id<-F.p.values(fit, term="time:oxygen")<0.001

##ASCA and APCA decomposition for every available term.
decomposition(asca, decomposition="pca", subset=id, scale="row")
decomposition(apca, decomposition="pca", subset=id, scale="row",
type="residual")

##Let's get the components for asca/apca decomposed objects
asca<-components(asca)
apca<-components(apca)

##Now let's try the PLSR decomposition for residuals and coefficients
plsr.residuals<-fit
plsr.coefficients<-fit
decomposition(plsr.coefficients, decomposition="plsr", subset=id,
scale="row")
decomposition(plsr.residuals, decomposition="plsr", subset=id, scale="row",
type="residual")

##Obtain the coefficients for decomposed plsr objects
##(coefficients/residuals)
plsr.coefficients<-components(plsr.coefficients)
plsr.residuals <- components(plsr.residuals)
}
}
\author{
  Cristobal Fresno and Elmer A Fernandez
}
\references{
  \enumerate{ \item Smilde AK, Jansen JJ, Hoefsloot HCJ,
  Lamer RAN, Van der Greef J, Timmerman ME (2005)
  ANOVA-simultaneous component analysis (ASCA): a new tool
  for analyzing designed meta-bolomics data, Bioinformatics
  21,13,3043 DOI:/10.1093/bioinformatics/bti476 \item
  Zwanenburg G, Hoefsloot HCJ, Westerhuis JA, Jansen JJ,
  Smilde AK (2011) ANOVA Principal component analysis and
  ANOVA-simultaneous component analysis: a comparison J.
  Chemometrics 25:561-567 DOI:10.1002/cem.1400 }
}
\seealso{
  \code{\link{prcomp}}, \code{\link{plsr}}
}

