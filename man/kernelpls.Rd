\name{kernelpls}
\alias{kernelpls}
\title{Kernel PLS (Dayal and MacGregor)}
\usage{kernelpls(X, Y, ncomp, newX)}
\description{This function should not be called directly, but through
  the generic \code{pls} function with the argument
  \code{method="kernel"} (default). Kernel PLS is particularly efficient 
  when the number of objects is (much) larger than the number of
  variables. The results are equal to the NIPALS algorithm. Several
  different forms of kernel PLS have been described in literature, e.g.
  by De Jong and Ter Braak, and two algorithms by Dayal and
  MacGregor. This function implements the
  fastest of the latter, not calculating the crossproduct matrix of
  X. In the Dyal & MacGregor paper, this is 'algorithm 1'.
  (\code{kernelpls}).
}
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
  \item{XvarExpl}{coeff. of determination (\code{X}).}
  \item{Ypred}{predictions for \code{newX} (if asked).}
}
\seealso{
  \code{\link{pls}}
  \code{\link{simpls}}
  \code{\link{mvr}}
}
\examples{
data(NIR)
attach(NIR)
NIR.kernelpls <- mvr(Xtrain, Ytrain, 1:6, validation="CV", method="kernelPLS")
}
\references{
  S. de Jong and C.J.F. ter Braak, J. Chemometrics, 8 (1994) 169-174
  B.S. Dayal and J. MacGregor, J. Chemometrics, 11 (1997) 73-85
  }
\keyword{regression}
