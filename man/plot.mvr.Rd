\name{plot.mvr}
\alias{plot.mvr}
\title{Plot a MVR object}
\usage{
plot.mvr(x,
         plottype=c("prediction", "validation", "coefficients",
                    "scores", "loadings"),
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
  \item{\code{type="loadings"}}{Shows a plot of loadings for X, in the
    case of a PCR or PLS1 model, and in the case of a PLS2
    model shows a biplot of loadings for X and Y. Loadings for Y are
    indicated in red.}
  \item{\code{type="scores"}}{Shows a plot of scores for X, in the
    case of a PCR model, and in the case of a PLS model shows a
    biplot of scores for X and Y. Scores for Y are
    indicated in red. In some cases, it may be necessary to scale scores
    before a meaningful plot is obtained.}
  \item{\code{type="coefficients"}}{Shows regression vector for the
    selected number of latent variables.}
  }
\arguments{
  \item{x}{an MVR object.} 
  \item{plottype}{type of plot. (Wow!)}
  \item{which}{whether training data (1), validation data (2) or both
    (1:2) should be shown in the \code{prediction} mode of this
    function.}
  \item{nlv}{the number of latent variables to take into account.  For the
    \code{loadings} and \code{scores} plotting types, this should be a
    vector of length 2, defaulting to c(1, 2). In the
    \code{coefficients} plotting type, this may be a vector. Not
    relevant for the \code{validation} plotting type.}
  \item{\dots}{several other arguments for the plotting command (not
    all!) can be used, as well. To see which ones, just try or look at
    the code...}
}
\section{References}{
  T. Hastie, R. Tibshirani and J. Friedman, \emph{The elements of
    statistical learning}, Springer (2001).
}
\examples{
data(NIR)
attach(NIR)
nir.pcr <- pcr(Xtrain,Ytrain,1:9, validat="CV")
\dontrun{plot(nir.pcr)}
plot(nir.pcr, "coefficients", nlv=3:6, type="l")
\dontrun{plot(nir.pcr, "prediction", nlv=3:6)}
plot(nir.pcr, "prediction", nlv=3:6, which=2)

data(sensory)
Pn <- scale(sensory$Panel)
Ql <- scale(sensory$Quality)
sens.simpls <- mvr(Ql, Pn, 1:5, method="SIMPLS")
plot(sens.simpls, "loadings")
plot(sens.simpls, "scores")
}
\keyword{regression}
