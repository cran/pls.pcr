\name{simpls}
\alias{simpls}
\title{De Jong's SIMPLS}
\usage{simpls(X, Y, ncomp, newX)}
\description{This function should not be called directly, but through
  the generic \code{pls} function with the argument
  \code{method="simpls"}. It is much faster than the NIPALS algorithm,
  especially when the number of X variables increases, but gives
  slightly different results in the case of multivariate Y. SIMPLS truly
  maximises the covariance criterion. According to De Jong, the standard
  PLS2 algorithms lie closer to ordinary least-squares regression where
  a precise fit is sought; SIMPLS lies closer to PCR with stable
  predictions.}
\arguments{
  \item{X}{a matrix of observations. \code{NA}s and \code{Inf}s are not
    allowed.} 
  \item{Y}{a vector or matrix of responses. \code{NA}s and \code{Inf}s
    are not allowed.} 
  \item{ncomp}{the number of latent variables to be used in the
    modelling. The default number of latent variables is the smallest of
    the number of objects or the number of variables in \code{X}.} 
  \item{newX}{optional new measurements: if present, predictions will be
    made for them.}
}
\value{A list containing the following
  components is returned:
  \item{B}{an array of regression coefficients for all items in
      \code{ncomp}. The dimensions of \code{B} are
      \code{c(nvar, npred, length(ncomp))} with \code{nvar} the number
      of \code{X} variables and \code{npred} the number of variables to be
      predicted in \code{Y}.}
    \item{XvarExpl}{Fraction of X-variance explained.}
    \item{YvarExpl}{Fraction of Y-variance explained (one column, even
      for multiple Y).
      }
  \item{Ypred}{predictions for \code{newX} (if asked).}
}
\seealso{
  \code{\link{kernelpls}}
  \code{\link{mvr}}
}
\examples{
data(NIR)
attach(NIR)
NIR.simpls <- mvr(Xtrain, Ytrain, 1:6, validation="CV", method="SIMPLS")
}
\references{S. de Jong, Chemometrics and Intelligent Laboratory
Systems, 18 (1993) 251-263.}
\keyword{regression}
