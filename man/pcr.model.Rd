\name{pcr.model}
\alias{pcr.model}
\title{PCR model building}
\usage{pcr.model(X, Y, ncomp, newX)}
\description{Function does the actual model building in PCR
  regression. Y is regressed on Principal Component scores.
  This function should not be called directly: the \code{pcr} 
  function is meant for that.}
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
\value{A list containing the following components is returned:
  \item{B}{an array of regression coefficients for all items in
      \code{ncomp}. The dimensions of \code{B} are
      \code{c(nvar, npred, length(ncomp))} with \code{nvar} the number
      of \code{X} variables and \code{npred} the number of variables to be
      predicted in \code{Y}.}
  \item{XvarExpl}{coeff. of determination (\code{X}).}
  \item{Ypred}{predictions for \code{newX} (if asked).}
}
\seealso{
  \code{\link{mvr}}
  \code{\link{pcr}}
}
\keyword{regression}
