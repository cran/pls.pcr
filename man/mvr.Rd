\name{mvr}
\alias{mvr}
\alias{pcr}
\alias{pls}
\title{Multivariate Regression}
\usage{
mvr(X, Y, ncomp=1:min(nobj, nvar),
    method=c("PCR", "SIMPLS", "kernelPLS"),
    validation=c("none","CV"), grpsize, niter)
pcr(\dots)
pls(\dots, method)
}
\description{Generic function for multivariate regression. PCR and two
  types of PLS (SIMPLS and kernel-PLS) are implemented. The functions
  \code{pcr} and \code{pls} are simple wrappers for \code{mvr}.
  Cross-validation can be used; a nmber of latent variables
  to retain is suggested.}
\arguments{
  \item{X}{a matrix of observations. \code{NA}s and \code{Inf}s are not
    allowed.} 
  \item{Y}{a vector or matrix of responses. \code{NA}s and \code{Inf}s
    are not allowed.} 
  \item{ncomp}{the numbers of latent variables to be assessed in the
    modelling. Default is from one to the rank of \code{X}.}
  \item{method}{the multivariate regression method to be used.}
  \item{validation}{validation method, either "none" or cross-validation
    ("CV").}
  \item{grpsize}{the group size for the "CV" validation. If not
    specified, this is determined by the \code{niter} parameter.}
  \item{niter}{the number of iterations in the
    cross-validation. Default: 10. If both \code{niter} and
    \code{grpsize} are specified, \code{niter} takes precedence.}
  \item{\dots}{arguments for \code{mvr}.}
}
\value{An object of class \code{mvr} is returned. This has the following
  components:
  \item{X}{original X data matrix.}
  \item{Y}{original Y data. Y is represented as a matrix, possibly with
    only one column.}
  \item{ncomp}{the number of latent variables considered in the PCR
    regression.}
  \item{method}{the method used.}
  \item{training}{a list with the following components:
    \item{B}{an array of regression coefficients for all items in
      \code{ncomp}. The dimensions of \code{B} are
      \code{c(nvar, npred, length(ncomp))} with \code{nvar} the number
      of \code{X} variables and \code{npred} the number of variables to be
      predicted in \code{Y}.
    }
    \item{Ypred}{predicted values for the training data.}
    \item{RMS}{Root-mean-square error of the training data for all items
      in \code{ncomp}.} 
    \item{R2}{coefficient of multiple determination for the training
      data.}
    \item{XvarExpl}{Fraction of X-variance explained.}
    \item{YvarExpl}{Fraction of Y-variance explained (one column, even
      for multiple Y).
    }
    \item{Xscores}{X-scores for objects.}
    \item{Xload}{X-loadings for objects.}
    \item{Yscores}{Y-scores for objects: only in case of a PLS model
      with more than one Y variable.}
    \item{Yload}{Y-loadings for objects: only in case of a PLS model
      with more than one Y variable.}
    }
    \item{validat}{a list with the following components:
      \item{Ypred}{predicted values.}
      \item{nLV}{number of latent variables suggested by
	cross-validation. The least complex model within one standard
	error of the best is chosen.
      }
      \item{niter}{the number of groups in the cross-validation. If this
	number equals the number of objects, this is leave-one-out
	validation.}
      \item{RMS}{Root-mean-square error from cross-validation for all items
	in \code{ncomp}.} 
      \item{RMS.sd}{standard error of the root-mean-square estimate.} 
      \item{R2}{coefficient of multiple determination from
	cross-validation.}
    }
}
\seealso{
  \code{\link{pcr.model}}, \code{\link{simpls}},
  \code{\link{kernelpls}}, \code{\link{plot.mvr}},
  \code{\link{summary.mvr}}}
}
\examples{
data(NIR)
attach(NIR)
NIR.pcr <- mvr(Xtrain, Ytrain, 1:6, validation="CV")
NIR.simpls <- mvr(Xtrain, Ytrain, 1:6, validation="CV", method="SIMPLS")
NIR.kernelpls <- mvr(Xtrain, Ytrain, 1:6, validation="CV", method="kernelPLS")

data(sensory)
Pn <- scale(sensory$Panel)
Ql <- scale(sensory$Quality)
sens.pcr <- mvr(Ql, Pn, 1:5)
sens.simpls <- mvr(Ql, Pn, 1:5, method="SIMPLS")
sens.kernelpls <- mvr(Ql, Pn, 1:5, method="kernelPLS")
}
\keyword{regression}
