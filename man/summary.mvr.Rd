\name{summary.mvr}
\alias{summary.mvr}
\title{Summary of a MVR object}
\usage{summary.mvr(object, what, digits=4, print.gap=3, ...)}
\description{The information in the MVR object is summarised: the
  function gives the dimensions of X and Y, the method used, and the
  number of latent variables considered. If \code{what} equals "all" or
  "training", the RMS and R2 values of the training set are given; if
  \code{what} equals "all" or "validation", RMS and R2 values for the
  cross-validation are given. In the latter case, a suggestion is made
  about the optimal number of latent variables, as suggested in Hastie
  et al. (2001).
}
\arguments{
  \item{object}{An MVR object.}
  \item{what}{one of "all", "training" or "validation". }
  \item{digits}{Number of digits, default is 4.}
  \item{print.gap}{gap between columns of the print table.}
  \item{\dots}{other arguments to \code{print.table}.}
}
\section{References}{
  T. Hastie, R. Tibshirani and J. Friedman, \emph{The elements of
    statistical learning}, Springer (2001).
}
\examples{
data(NIR)
attach(NIR)
nir.mvr <- mvr(Xtrain, Ytrain, 2:8, valid="CV", grps=1)
summary(nir.mvr)
}
\keyword{regression}
