plot.mvr <- function(x,
                     plottype=c("prediction", "validation", "coefficients"),
                     nlv=mvrmodel$validat$nLV,
                     which=1:2, ...)
{
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  
  mvrmodel <- x
  npred <- dim(mvrmodel$Y)[2]
  if (is.null(ynames <- dimnames(mvrmodel$Y)[[2]]))
    ynames <- ifelse (npred > 1, list(paste("Y", 1:npred)), list("Y"))[[1]]
  
  plottype <- match.arg(plottype)
  if (is.null(mvrmodel$validat)) {
    if (plottype == "validation")
      stop("No validation data available")
    if (plottype == "prediction")
      which = 1
  }

  if (is.null(nlv))
    nlv <- mvrmodel$ncomp
  
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
                       ifelse (mvrmodel$validat$niter == dim(mvrmodel$X)[1],
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
           if (length(nlv) > 1) {
             nrow <- floor(sqrt(length(nlv)))
             mfrow <- c(nrow, ceiling(length(nlv)/nrow))
             par(mfrow=mfrow)
           }
           if (npred > 1) par(ask = TRUE)
           
           for (i in 1:npred) {
             for (j in seq(along=nlv)) {
               plot(mvrmodel$training$B[ , i, nlv[j]], xlab="Variable",
                    ylab="Regression coefficient",
                    main=paste("Regression vector using", nlv[j],
                      "latent variables"), ...)
               if (npred > 1) mtext(ynames[i])
             }
           }
         },
         prediction = {
           show <- rep(FALSE, 2)
           if (!is.numeric(which) || any(which < 1) || any(which > 2))
             stop("`which' must be in 1:2")
           show[which] <- TRUE
           if (length(which)>1) {
             par(mfrow = c(length(nlv),2))
           } else {
             if (length(nlv)> 1) {
               nrow <- floor(sqrt(length(nlv)))
               mfrow <- c(nrow, ceiling(length(nlv)/nrow))
               par(mfrow=mfrow)
             }
           }
          
           if (npred > 1)
             par(ask = TRUE )
           par(pty="s")
           
           for (i in 1:npred) {
             for (j in seq(along=nlv)) {
               if (show[1]) {
                 plot(mvrmodel$Y[,i], mvrmodel$training$Ypred[,i,nlv[j]],
                      xlab = "Measured", ylab = "Predicted", type = "n",
                      xlim = range(mvrmodel$Y[,i],
                        mvrmodel$training$Ypred[,i,nlv[j]]),
                      ylim = range(mvrmodel$Y[,i],
                        mvrmodel$training$Ypred[,i,nlv[j]]),
                      main=paste("Training data\n",nlv[j],"latent variables"))
                 if (npred > 1) mtext(ynames[i])
                 abline(0, 1, col="blue")
                 points(mvrmodel$Y[,i], mvrmodel$training$Ypred[,i,nlv[j]])
               }
               
               if (show[2]) {
                 plot(mvrmodel$Y[,i], mvrmodel$validat$Ypred[,i,nlv[j]],
                      xlab="Measured", ylab="Predicted", type="n",
                      xlim=range(mvrmodel$Y[,i],
                        mvrmodel$validat$Ypred[,i,nlv[j]]),
                      ylim=range(mvrmodel$Y[,i],
                        mvrmodel$validat$Ypred[,i,nlv[j]]),
                      main=paste("Cross-validation data\n",
                        nlv[j],"latent variables"))
                 if (npred > 1) mtext(ynames[i])
                 abline(0, 1, col="blue")
                 points(mvrmodel$Y[,i], mvrmodel$validat$Ypred[,i,nlv[j]])
               }
             }
           }
         }
         )
}
