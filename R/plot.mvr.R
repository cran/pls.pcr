plot.mvr <- function(x,
                     plottype=c("prediction", "validation",
                       "coefficients", "scores", "loadings"),
                     nlv=NULL, which=1:2, ...)
{
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  mvrmodel <- x
  npred <- mvrmodel$npred
  ynames <- dimnames(mvrmodel$Y)[[2]]
  
  plottype <- match.arg(plottype)
  if (is.null(mvrmodel$validat)) {
    if (plottype == "validation")
      stop("No validation data available")
    if (plottype == "prediction")
      which = 1
  }

  switch(plottype,
         validation = {
           if (npred > 1) par(ask = TRUE)
           
           for (i in 1:npred) {
             x <- mvrmodel$ncomp
             y <- matrix(c(mvrmodel$validat$RMS[,i],
                           mvrmodel$validat$RMS[,i]
                           + mvrmodel$validat$RMS.sd[,i],
                           mvrmodel$validat$RMS[,i]
                           - mvrmodel$validat$RMS.sd[,i]),
                         ncol=3)
             
             matplot(x, y, xlab="# LV's", ylab="CV error",
                     main=paste(mvrmodel$method, ": ",
                       ifelse (mvrmodel$validat$niter == mvrmodel$nobj,
                               "LOO-CV",
                               list(paste(mvrmodel$validat$niter, "-fold CV",
                                          sep=""))[[1]]),
                       sep=""),
                     type="n", ...)
             if (npred > 1) mtext(ynames[i])
             abline(h=min(y[,2]), lty=3)
             abline(v=mvrmodel$validat$nLV, lty=3)
             lines(x, y[,1], col=2)
             segments(x, y[,2], x, y[,3], col=3)
             segments(c(x-.1, x-.1), c(y[,2], y[,3]),
                      c(x+.1, x+.1), c(y[,2], y[,3]), col=3)
             points(x, y[,1], col=2, pch=16, cex=1.2)
           }
         },
         coefficients = {
           if (is.null(nlv))
             nlv <- mvrmodel$valid$nLV
           if (is.null(nlv))
             nlv <- mvrmodel$ncomp
  
           if (length(nlv) > 1) {
             nrow <- floor(sqrt(length(nlv)))
             mfrow <- c(nrow, ceiling(length(nlv)/nrow))
             par(mfrow=mfrow)
           }
           if (npred > 1) par(ask = TRUE)
           
           for (i in 1:npred) {
             for (j in seq(along=nlv)) {
               plot(mvrmodel$training$B[ , i,
                                        which(mvrmodel$ncomp == nlv[j])],
                    xlab="Variable", ylab="Regression coefficient",
                    main=paste("Regression vector using", nlv[j],
                      "latent variables"), ...)
               abline(h=0, col="gray")
               if (npred > 1) mtext(ynames[i])
             }
           }
         },
###         inner = {
###           if (length(nlv) > 1) {
###             nrow <- floor(sqrt(length(nlv)))
###             mfrow <- c(nrow, ceiling(length(nlv)/nrow))
###             par(mfrow=mfrow)
###           }
###           if (npred > 1) par(ask = TRUE)
###         },
         scores = {
           if (is.null(nlv) | length(nlv) != 2) nlv = 1:2

           if (mvrmodel$method != "PCR") {
             biplot(mvrmodel$training$Xscores[,nlv],
                    mvrmodel$training$Yscores[,nlv],
                    var.axes=FALSE, main="Scores",
                    xlab=paste("LV", nlv[1]),
                    ylab=paste("LV", nlv[2]), ...)
             abline(h=0, v=0, col="gray")
           } else { # only X scores
             plot(mvrmodel$training$Xscores[,nlv],
                  main="Scores",
                  xlab=paste("LV", nlv[1]),
                  ylab=paste("LV", nlv[2]), type="n", ...)
             abline(h=0, v=0, col="gray")
             labs <- dimnames(mvrmodel$Xscores)[[1]]
             if (is.null(labs)) labs <- 1:mvrmodel$nobj
             text(mvrmodel$training$Xscores[,nlv[1]],
                  mvrmodel$training$Xscores[,nlv[2]],
                  labs)
           }
         },
         loadings = {
           if (is.null(nlv) | length(nlv) != 2) nlv = 1:2
           if (mvrmodel$method != "PCR" && mvrmodel$npred > 1) {
             biplot(mvrmodel$training$Xload[,nlv],
                    mvrmodel$training$Yload[,nlv],
                    var.axes=FALSE, main="Loadings",
                    xlab=paste("LV", nlv[1]),
                    ylab=paste("LV", nlv[2]), ...)
             abline(h=0, v=0, col="gray")
           } else { # only X loadings
             plot(mvrmodel$training$Xload[,nlv],
                  main="Loadings",
                  xlab=paste("LV", nlv[1]),
                  ylab=paste("LV", nlv[2]), type="n", ...)
             abline(h=0, v=0, col="gray")
             labs <- dimnames(mvrmodel$Xload)[[1]]
             if (is.null(labs)) labs <- 1:mvrmodel$nobj
             text(mvrmodel$training$Xload[,nlv[1]],
                  mvrmodel$training$Xload[,nlv[2]],
                  labs)
           }
         },
         prediction = {
           if (is.null(nlv))
             nlv <- mvrmodel$valid$nLV
           if (is.null(nlv))
             nlv <- mvrmodel$ncomp
  
           show <- rep(FALSE, 2)
           ### which == 1: show predictions of training set
           ### which == 2: show crossvalidation predictions
           ### which == c(1,2): show both
           if (!is.numeric(which) || any(which < 1) || any(which > 2))
             stop("`which' must be in 1:2")
           show[which] <- TRUE
           if (length(which)>1) {
             ## show at most three rows of two plots
             ## if both training and validation plots are required
             if (length(nlv) > 3) {
               par("ask" = TRUE)
               par(mfrow = c(3,2))
               nScreens <- ceiling(length(nlv) / 3)
             } else {
               nScreens <- 1
               par(mfrow = c(length(nlv),2))
             }
           } else {
             ## if only validation or training plots, then cram
             ## everything into one screen
             if (length(nlv)> 1) {
               nrow <- floor(sqrt(length(nlv)))
               mfrow <- c(nrow, ceiling(length(nlv)/nrow))
               par(mfrow=mfrow)
             }
             nScreens <- 1
           }
         
           if (npred > 1)
             par(ask = TRUE )
           par(pty="s")
           
           for (i in 1:npred) {
             allindices <- match(nlv, mvrmodel$ncomp)
             if (any(is.na(allindices)))
               stop(paste("Incorrect number of LVs; should be a subset of",
                          paste(mvrmodel$ncomp, collapse=", ")))
             for (k in 1:nScreens) {
               if (nScreens == 1) {
                 indices <- allindices
               } else {
                 cnts <- (k-1)*3 + 1:3
                 cnts <- cnts[cnts <= length(nlv)]
                 indices <- allindices[cnts]
               }
               for (j in indices) {
                 if (show[1]) {
                   plot(mvrmodel$Y[,i], mvrmodel$training$Ypred[,i,j],
                        xlab = "Measured", ylab = "Predicted", type = "n",
                        xlim = range(mvrmodel$Y[,i],
                          mvrmodel$training$Ypred[,i,j]),
                        ylim = range(mvrmodel$Y[,i],
                          mvrmodel$training$Ypred[,i,j]),
                        main=paste("Training data\n",
                          dimnames(mvrmodel$training$Ypred)[[3]][j]))
                   if (npred > 1) mtext(ynames[i])
                   abline(0, 1, col="blue")
                   points(mvrmodel$Y[,i], mvrmodel$training$Ypred[,i,j])
                 }
                 
                 if (show[2]) {
                   plot(mvrmodel$Y[,i], mvrmodel$validat$Ypred[,i,j],
                        xlab="Measured", ylab="Predicted", type="n",
                        xlim=range(mvrmodel$Y[,i],
                          mvrmodel$validat$Ypred[,i,j]),
                        ylim=range(mvrmodel$Y[,i],
                          mvrmodel$validat$Ypred[,i,j]),
                        main=paste("Cross-validation data\n",
                          dimnames(mvrmodel$validat$Ypred)[[3]][j]))
                   if (npred > 1) mtext(ynames[i])
                   abline(0, 1, col="blue")
                   points(mvrmodel$Y[,i], mvrmodel$validat$Ypred[,i,j])
                 }
               }
             }
           }
         }
         )
}
