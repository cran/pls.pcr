summary.mvr <- function(object, what=c("all", "validation", "training"),
                        digits=4, print.gap=3, ...)
{
  what <- match.arg(what)
  if (what == "all") what <- c("validation", "training")
  if (is.null(object$validat)) what <- "training"
  
  cat("Data: \tX dimension:", object$nobj, object$nvar,
      "\n\tY dimension:", object$nobj, object$npred)
  cat("\nMethod:", object$method)
  cat("\nNumber of latent variables considered:",
      ifelse(length(object$ncomp) > 1 && max(diff(object$ncomp) == 1),
             paste(min(object$ncomp), "-", max(object$ncomp), sep=""),
             object$ncomp))
  if ("validation" %in% what)
    cat("\nSuggested number of latent variables:", object$validat$nLV,
        "(cross validation)\n")
  
  for (wh in what) {
    if (wh == "training") {
      obj <- object$training
      cat("\n\nTRAINING: % variance explained\n")
    } else {
      obj <- object$validat
      cat("\n\nVALIDATION:\n")
    }

    if (wh == "training") {
      xve <- diag(crossprod(obj$Xload)) / object$Xss
      yve <- matrix(0, length(object$ncomp), object$npred)
      for (i in seq(along=object$ncomp))
        for (j in 1:object$npred)
          yve[i,j] <- cor(object$Y[,j], obj$Ypred[,j,i])^2
      
      tbl <- cbind(100*cumsum(xve)[object$ncomp], 100*yve)
      dimnames(tbl) <- list(dimnames(obj$RMS)[[1]],
                            c("X", dimnames(object$Y)[[2]]))
      print(tbl, format="f", print.gap=print.gap, digits=digits, ...)
    } else {
      for (i in 1:object$npred) {
        cat(dimnames(obj$RMS)[[2]][i],"\n")
        
        tbl <- matrix(0, length(object$ncomp), 3)
        dimnames(tbl) <- list(dimnames(obj$RMS)[[1]],
                              c("RMS", "sd(RMS)", "Q^2"))
        tbl[,1] <- obj$RMS[,i]
        tbl[,2] <- obj$RMS.sd[,i]
        tbl[,3] <- obj$R2[,i]

        print(tbl, format="f", print.gap=print.gap, digits=digits, ...)
        cat("\n")
      }
    }
  }
}

