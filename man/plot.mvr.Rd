\name{plot.mvr}
\alias{plot.mvr}
\title{Plot a MVR object}
\usage{
plot.mvr(x,
         plottype=c("validation", "prediction", "coefficients"),
         nlv, which=1:2, \dots)
}
\description{Several different ways of plotting a \code{mvr} object are
  supported:
  \item{\code{type="validation"}}{plot of errors versus the number of
    latent variables. Standard errors are also shown, and the selected
    number of components is indicated, as suggested in Hastie et
    al. (2001).}
  \item{\code{type="prediction"}}{For each Y variable, a plot is
    produced of true versus predicted values. If validation data are
    available, both training and validation plots are shown.
    For each Y variable separate plots are made.
  }
  \item{\code{type="coefficients"}}{Shows regression vector for the
    selected number of latent variables.}
  }
\arguments{
  \item{x}{an MVR object.} 
  \item{plottype}{Type of plot. (Wow!)}
  \item{which}{Whether training data (1), validation data (2) or both
    (1:2) should be shown in the \code{prediction} mode of this
    function.}
  \item{nlv}{the number of latent variables (in the \code{coefficients}
    plotting type, this may be a vector) to be used for the
    plot. Not relevant for the \code{validation} plotting type.} 
  \item{\dots}{optional arguments for the plotting command.}
}
\section{References}{
  T. Hastie, R. Tibshirani and J. Friedman, \emph{The elements of
    statistical learning}, Springer (2001).
}
\examples{
data(NIR)
attach(NIR)
nir.pcr <- pcr(Xtrain,Ytrain,1:9, validat="CV")
plot(nir.pcr)
plot(nir.pcr, "coefficients", nlv=3:6, type="l")
plot(nir.pcr, "prediction", nlv=3:6)
plot(nir.pcr, "prediction", nlv=3:6, which=2)
}
\keyword{regression}
