\name{predict.mvr}
\alias{predict.mvr}
\title{MVR predictions}
\usage{predict.mvr(object, newX, nlv, \dots)}
\description{Prediction for MVR (PCR, PLS) models. New responses are
  predicted using a previously defined model and a new matrix of
  independent observations.}
\arguments{
  \item{object}{a MVR model} 
  \item{newX}{a matrix of new independent observations} 
  \item{nlv}{the number of latent variables (optional)}
  \item{\dots}{not used currently.}
}
\value{
  Predicted values are returned.
}
\seealso{
  \code{\link{plot.mvr}}
  \code{\link{summary.mvr}}
}
\examples{
data(NIR)
attach(NIR)
nir.mvr <- mvr(Xtrain, Ytrain, 2)
nir.predictions <- predict(nir.mvr, Xtest)
}
\keyword{regression}
