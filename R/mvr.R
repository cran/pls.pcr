cat("\nThis package has been superseded by the package",
    "\n'pls', authored by Bjørn-Helge Mevik and Ron Wehrens.",
    "\nPackage 'pls.pcr' will not be actively supported anymore,",
    "\nand will eventually be deprecated.\n\n")


### Generic function for PCR and PLS.
### Contains preprocessing and crossvalidation material; real
### modelling is done in separate functions.

mvr <- function(X, Y, ncomp,
                method=c("PCR", "SIMPLS", "kernelPLS"), 
                validation=c("none","CV"),
	        grpsize, niter) 
{
  if (is.vector(Y))
    Y <- matrix(Y, ncol=1)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  nobj <- dim(X)[1]
  nvar <- dim(X)[2]
  npred <- dim(Y)[2]

  objnames <- dimnames(X)[[1]]
  if (is.null(objnames)) objnames <- dimnames(Y)[[1]]
  ### if too many objects, names are a nuisance so dont use the next line
  ##  if (is.null(objnames)) objnames <- paste("obj", 1:nobj, sep=".")
  xvarnames <- dimnames(X)[[2]]
  ### if too many variables, names are a nuisance too
  ##  if (is.null(xvarnames)) xvarnames <- paste("X", 1:nvar, sep=".")
  yvarnames <- dimnames(Y)[[2]]
  ## not often that you have too many Y variables...
  if (is.null(yvarnames))
    yvarnames <- ifelse (npred > 1,
                         paste("Y", 1:npred, sep="."),
                         "Y")
  dimnames(X) <- dimnames(Y) <- NULL

  validation <- match.arg(validation)
  if (validation == "CV") {
    if (missing(grpsize) & missing(niter)) {
      if (nobj > 20) {
        niter <- 10
        grpsize <- nobj %/% niter
      } else {
        niter <- nobj
        grpsize <- 1
      }
    } else {
      if (missing(niter))
        niter <- nobj %/% grpsize
      else
        grpsize <- nobj %/% niter
    }
    if (niter*grpsize == nobj)
      maxlv <- min(nobj-grpsize, nvar)
    else
      maxlv <- min(nobj-grpsize-1, nvar)
  } else {
    maxlv <- min(nobj, nvar)
  }
  
  method <- match.arg(method)
  modelfun <- switch(method,
                     PCR = pcr.model,
                     SIMPLS = simpls,
                     kernelPLS = kernelpls)

  if (missing(ncomp))
    ncomp <- 1:maxlv

  if (max(ncomp) > maxlv) {
    ncomp <- ncomp[ncomp<=maxlv]
    if (is.null(ncomp))
      stop("Invalid number of latent variables, ncomp")
    else
      warning("Reset maximum number of latent variables to ", max(ncomp))
  }
  ncompnames <- paste(ncomp, "LV's")
  
  Xm <- scale(X, scale=FALSE)
  Ym <- scale(Y, scale=FALSE)
  training <- c(modelfun(Xm, Ym, ncomp, Xm),
                list(RMS = matrix(0, length(ncomp), npred)))
    
  for (i in seq(along=ncomp)) {
    training$Ypred[ , , i] <-
      sweep(training$Ypred[ , , i, drop=FALSE], 2,
            colMeans(Y), FUN = '+')
    training$RMS[i,] <- sqrt(colMeans((Y-training$Ypred[ , , i])^2))
  }

  ## Not ncompnames here, since all scores and loadings from 1:max(ncomp)
  ## are calculated and stored, not only for those numbers in vector ncomp.
  dimnames(training$Xscores) <- list(objnames,
                                     paste("LV", 1:max(ncomp)))
  dimnames(training$Xload) <- list(xvarnames,
                                   paste("LV", 1:max(ncomp)))
  if (method != "PCR") {
    dimnames(training$Yscores) <- list(objnames,
                                       paste("LV", 1:max(ncomp)))
    dimnames(training$Yload) <- list(yvarnames,
                                     paste("LV", 1:max(ncomp)))
  }
  dimnames(training$RMS) <- list(ncompnames, yvarnames)
  dimnames(training$Ypred) <- list(objnames, yvarnames, ncompnames)
  dimnames(Y) <- list(objnames, yvarnames)
  
  mvrmodel <- list(nobj=nobj, nvar=nvar, npred=npred, ncomp=ncomp,
                   Y = Y, Xmeans=attr(Xm, "scaled:center"),
                   Xss = ifelse(nobj<nvar,
                     sum(apply(Xm, 1, crossprod)),
                     sum(apply(Xm, 2, crossprod))),
                   training=training, method=method)
  
  if (validation == "CV") {
    validat <- list(niter = niter, nLV=-1,
                    RMS = matrix(0, length(ncomp), npred),
                    RMS.sd = matrix(0, length(ncomp), npred),
                    R2 = matrix(0, length(ncomp), npred))
    
    Ypred <- array(0, c(nobj, npred, length(ncomp)))
    indices <- sample(1:nobj)
    setsizes <- rep(grpsize, niter)
    if (grpsize*niter < nobj)
      setsizes[1:(nobj-niter*grpsize)] <- grpsize + 1
    lastone <- 0
    errors <- matrix(0, nobj, ncol(Y))
    for (i in 1:niter) {
      start <- lastone + 1
      end <- lastone + setsizes[i]
      xmns <- colMeans(X[-indices[start:end],])
      ymns <- colMeans(Y[-indices[start:end],,drop=FALSE])
      Xv <- sweep(X[-indices[start:end],], 2, xmns)
      Xv.test <- sweep(X[indices[start:end],,drop=FALSE], 2, xmns)
      Yv <- sweep(Y[-indices[start:end],,drop=FALSE], 2, ymns)
      
      Ypred[indices[start:end], , ] <-
        sweep(modelfun(Xv, Yv, ncomp, Xv.test)$Ypred, 2, ymns, FUN='+')
      
      lastone <- end
    }
    
    for (i in seq(along=ncomp)) {
      Ydiff <- Y - Ypred[,,i]
      validat$RMS[i,] <- sqrt(colMeans(Ydiff^2))
      validat$RMS.sd[i,] <- apply(Ydiff, 2, sd)/sqrt(nobj)
      validat$R2[i,] <- diag(cor(Y, Ypred[,,i]))^2
    }
    if (npred > 1) {
      rmsmat <- rowSums(validat$RMS)
      rmsmat.sd <- rowSums(validat$RMS.sd)
      validat$nLV <- ncomp[min(which(rmsmat < min(rmsmat+rmsmat.sd)))]
    } else {
      validat$nLV <-
        ncomp[min(which(validat$RMS < min(validat$RMS + validat$RMS.sd)))]
    }
    
    dimnames(validat$RMS) <- dimnames(validat$RMS.sd) <-
      dimnames(validat$R2) <- list(ncompnames, yvarnames)
    dimnames(Ypred) <- dimnames(mvrmodel$training$Ypred)
    validat$Ypred <- Ypred
    mvrmodel <- c(mvrmodel, list(validat=validat))
  }

  class(mvrmodel) <- "mvr"
  mvrmodel
}
