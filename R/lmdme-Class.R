#' lmdme S4 class: Linear Model decomposition for Designed Multivariate
#' Experiments. 
#'
#' Linear Model ANOVA decomposition of Designed Multivariate Experiments based
#' on limma \code{\link{lmFit}} implementation. For example in a two factor
#' experimental design with interaction, the linear model of the i-th
#' observation (gene) can be written: 
#' \cr \eqn{X=\mu+A+B+AB+\epsilon} \cr
#' where 
#' \itemize{ 
#'  \item X stands for the observed value 
#'  \item The intercept \eqn{\mu} 
#'  \item A, B and AB are the first, second and interaction terms respectively
#'  \item The error term \eqn{\epsilon ~ N(0,\sigma^2)}.
#'  } 
#' The the model is iterative decomposed in a step by step fashion decomposing
#' one term at each time: 
#' \enumerate{ 
#'  \item The intercept is estimated using \eqn{X=\mu+E_1} 
#'  \item The first factor (A) using \eqn{E_1=A+E_2} 
#'  \item The second factor (B) using \eqn{E_2=B+E_3} 
#'  \item The interaction (AB) using \eqn{E_3=AB+E_4}.
#' } 
#' For each decomposed step the model, residuals, coefficients, p-values and
#' F-values are stored in a list container, so their corresponding length is
#' equal to the number of model terms + 1 (the intercept). 
#'
#' @section Features: 
#' \enumerate{
#'   \item Flexible formula type interface, 
#'   \item Fast limma based implementation based on \code{\link{lmFit}},
#'   \item p values for each estimated coefficient levels in each factor
#'   \item F values for factor effects 
#'   \item Plotting functions for PCA and PLS.
#' }
#'
#' @section Slots: 
#' \itemize{
#'  \item design: data.frame with experimental design.
#'  \item model: formula with the designed model to be decomposed.
#'  \item modelDecomposition: list with the model formula obtained for each 
#'  deflation step.
#'  \item residuals: list of residual matrices G rows(genes) x N columns
#'  (arrays-designed measurements).
#'  \item coefficients: list of coefficient matrices. Each matrix will have 
#'  G rows(genes) x k columns(levels of the corresponding model term).
#'  \item p.values: list of p-value matrices.
#'  \item F.p.values: list with corresponding F-p-values vectors for each 
#'  individual.
#'  \item components: list with corresponding PCA or PLS components for the 
#'  selected term/s.
#'  \item componentsType: name character vector to keep process trace of the 
#'  variance/covariance components slot: decomposition ("pca" or "pls"), type
#'  ("apca" for ANOVA-PCA or "asca" for ANOVA-SCA) and scale ("none", "row" or
#'  "column") 
#' }
#'
#' @section lmdme-general-functions:
#' \describe{
#'  \item{print, show}{Basic output for lmdme class}
#'  \item{summary}{Basic statistics for lmdme class} 
#'  \item{design, model, modelDecomposition, residuals and coefficients}{Getters
#' for their respective slots.} 
#'  \item{p.values, F.p.values, components and componentsType}{Getters for their
#' respective slots.}
#' }
#'
#' @section ANOVA-linear-decomposition-functions:
#' \describe{
#'  \item{lmdme}{Function that produce the complete ANOVA decomposition based 
#'  on model specification through formula interface. Technically is a high
#'  level wrapper of initialize function.}
#'  \item{modelDecomposition}{Getter for decomposed used formula in each step}
#'  \item{p.adjust}{Adjust coefficients p-values for Multiple Comparisons Test.}
#'  \item{Fpvalues, pvalues}{Getters for corresponding model decomposed 
#'  associated coefficient statistics in each step, for each observation}
#'  \item{residuals, resid, coef, coefficients, fitted.values, fitted}{Getters
#'   for corresponding model decomposed in each step}
#'  \item{permutation}{Produces de specified lmdme in addition to the required
#'  permuted objects (sampling the columns of data), using the same parameters
#'  to fit the model.}
#' }
#'
#' @section variance-covariance-decomposition-functions:
#' \describe{
#'  \item{decomposition}{Function to perform PCA or PLS over the ANOVA
#'  decomposed terms. PCA can be performed over \eqn{E_1}, \eqn{E_2} or
#'  \eqn{E_3} and it is referred as ANOVA-PCA (APCA) but, if it is
#'  performed over the coefficients it is referred as ANOVA-SCA (ASCA). On
#'  the other hand PLSR is based on pls library and if it is performed on
#'  coefficients (ASCA like) it uses the identity matrix for output
#'  co-variance maximization or can be carried out over the \eqn{E_{1,2 or 3}}
#'  (APCA like) using the design matrix as output.}
#'  \item{components}{Getter for PCA or PLS decomposed models.}
#'  \item{componentsType}{Getter for componentsType slot.}
#'  \item{leverage}{Leverage calculation over PCA (APCA or ASCA) terms.}
#'  \item{biplot}{Biplots for PCA or PLSR decomposed terms.}
#'  \item{screeplot}{Screeplot over each decomposed PCA decomposed term.}
#'  \item{loadingplot}{Loadingplot for PCA interaction terms.}
#' }
#'
#' @author Cristobal Fresno and Elmer A Fernandez
#' 
#' @seealso \code{\link{lmdme}}, \code{\link{decomposition}}, 
#'  \code{\link{biplot}}, \code{\link{loadingplot}} and additional related
#'  lmdme class functions.
#'
#' @references 
#' \enumerate{
#'  \item Smilde AK, Jansen JJ, Hoefsloot HCJ, Lamer RAN, Van der Greef J, 
#'  Timmerman ME (2005) ANOVA-simultaneaus component analysis (ASCA): a new tool
#'  for analyzing designed \cr metabolomics data, Bioinformatics 21,13,3043
#'  DOI:/10.1093/bioinformatics/bti476
#'  \item Zwanenburg G, Hoefsloot HCJ, Westerhuis JA, Jansen JJ, Smilde AK
#'  (2011) ANOVA.Principal component analysis and ANOVA-simultaneaus component
#'  analysis: a comparison J. \cr Chemometrics 25:561-567 DOI:10.1002/cem.1400
#'  \item Tarazona S, Prado-Lopez S, Dopazo J, Ferrer A, Conesa A (2012) 
#'  Variable Selection for Multifactorial Genomic Data, Chemometrics and
#'  Intelligent Laboratory Systems, 110:113-122
#' }
#'
#' @name lmdme-class
#' @rdname lmdme-Class
#' @exportClass lmdme
setClass(Class="lmdme",
  representation=representation(
    design="data.frame",
    model="formula",
    modelDecomposition="list",
    residuals="list",
    coefficients="list",
    p.values="list",
    F.p.values="list",
    components="list",
    componentsType="character"),
  prototype=prototype(
    design=data.frame(),
    residuals=list(),
    coefficients=list(),
    p.values=list(),
    F.p.values=list(),
    modelDecomposition=list(),
    components=list(),
    componentsType=character()),
  validity=function(object){
    ## Check length of list of residuals, coefficients, pvalues, Fpvalues, models
    modelsLength<-c(length(object@residuals), length(object@coefficients), 
      length(object@p.values), length(object@F.p.values),
      length(object@modelDecomposition))
    if(any(modelsLength != modelsLength[1])) 
      stop("List type slots do not have the appropriate length")
    return(TRUE)
  }
)

