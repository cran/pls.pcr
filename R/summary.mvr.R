summary.mvr <- function(object, what=c("all", "training", "validation"), ...)
{
  what <- match.arg(what)
  if (is.null(object$validat))
    what <- "training"
  
  cat("Data: \tX dimension:", dim(object$X),
      "\n\tY dimension:", dim(object$Y))
  cat("\nMethod:", object$method)
  cat("\nNumber of latent variables considered:",
      ifelse(length(object$ncomp) > 1 & max(diff(object$ncomp) == 1),
             paste(min(object$ncomp), "-", max(object$ncomp), sep=""),
             object$ncomp), "\n")

  if (what == "all" | what == "training") {
    cat("\n\nTRAINING:\nRMS table:\n")
    print(object$training$RMS, format="f", digits=3)
    
    cat("\nCumulative fraction of variance explained:\n")
    varexpl <- cbind(object$training$XvarExpl, object$training$YvarExpl)
    print(varexpl, format="f", digits=3)
  }
  
  if (what == "all" | what == "validation") {
    cat("\n\nVALIDATION\nOptimal number of latent variables:",
        object$validat$nLV)
    cat("\n\nRMS table (",
        object$validat$niter, "-fold crossvalidation):\n", sep="")
    print(object$validat$RMS, format="f", digits=3)
    cat("\nCoefficient of multiple determination (R2):\n")
    print(object$validat$R2, format="f", digits=2)
    cat("\n")
  }
}
