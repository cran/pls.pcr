summary.mvr <- function(object, what=c("all", "validation", "training"), ...)
{
  what <- match.arg(what)
  if (what == "all") what <- c("validation", "training")
  if (is.null(object$validat)) what <- "training"
  
  cat("Data: \tX dimension:", dim(object$X),
      "\n\tY dimension:", dim(object$Y))
  cat("\nMethod:", object$method)
  cat("\nNumber of latent variables considered:",
      ifelse(length(object$ncomp) > 1 & max(diff(object$ncomp) == 1),
             paste(min(object$ncomp), "-", max(object$ncomp), sep=""),
             object$ncomp))
  if ("validation" %in% what)
    cat("\nSuggested number of latent variables:", object$validat$nLV,
        "(cross validation)\n")
  
  for (wh in what) {
    if (wh == "training") {
      obj <- object$training
      cat("\n\nTRAINING:\n")
    } else {
      obj <- object$validat
      cat("\n\nVALIDATION:\n")
    }

    if (wh == "training") {
      cat("Variance explained (%):\n")
      tbl <- cbind(obj$XvarExpl*100, obj$YvarExpl*100)
      print(tbl, format="f", print.gap=3, ...)
      cat("\nRMS:\n")
      print(obj$RMS, print.gap=3, ...)
    } else {
      for (i in 1:ncol(object$Y)) {
        cat(dimnames(obj$RMS)[[2]][i],"\n")
        
        tbl <- matrix(0, length(object$ncomp), 3)
        dimnames(tbl) <- list(dimnames(obj$RMS)[[1]],
                              c("RMS", "sd(RMS)", "Q^2"))
        tbl[,1] <- obj$RMS[,i]
        tbl[,2] <- obj$RMS.sd[,i]
        tbl[,3] <- obj$R2[,i]

        print(tbl, format="f", print.gap=3, ...)
        cat("\n")
      }
    }
  }
}

