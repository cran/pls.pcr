predict.mvr <- function(object, newX, nlv, ...)
{
  X <- sweep(newX, 2, object$Xmeans)

  if (length(object$ncomp) == 1) {
    if (!missing(nlv))
      if (nlv != object$ncomp)
        warning("overriding nlv argument: set to", object$ncomp)
    nlv <- object$ncomp
  } else {
    if (missing(nlv)) {
      if (is.null(object$validat))
        stop("no number of latent variables specified")
      else
        nlv <- object$validat$nLV
    }
  }

  index <- which(object$ncomp == nlv)
  sweep(X %*% object$training$B[ , , index, drop=TRUE], 2,
        colMeans(object$Y), FUN="+")
}
