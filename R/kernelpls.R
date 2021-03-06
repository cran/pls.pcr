### Contains R/Splus code for Kernel PLS 
### The method is particularly efficient if the number of
### variables is much smaller than the number of objects.
### Refined method by Dayal and MacGregor, J. Chemometr. 11 73-85, (1997).
### Modified kernel algorithm 1

kernelpls <- function(X, Y, ncomp, newX)
{
  nobj <- dim(X)[1]
  nvar <- dim(X)[2]
  npred <- dim(Y)[2]
  
  PP <- matrix(0, ncol=max(ncomp), nrow=nvar) # X loadings
  QQ <- matrix(0, ncol=max(ncomp), nrow=npred)# Y loadings
  TT <- matrix(0, ncol=max(ncomp), nrow=nobj)
  RR <- matrix(0, ncol=max(ncomp), nrow=nvar)
  UU <- matrix(0, ncol=max(ncomp), nrow=nobj)
  B <- array(0, c(dim(X)[2], dim(Y)[2], length(ncomp)))

  if (!is.null(newX))
    Ypred <- array(0, c(dim(newX)[1], npred, length(ncomp)))
                        
  XtY <- crossprod(X, Y)
  for (a in 1:max(ncomp)) {
    if (npred==1)
      ww <- XtY
    else {
      if (npred>nvar)
        qq <- La.svd(XtY)$vt[1,]
      else
        qq <- La.svd(XtY)$u[1,]
      ww <- XtY %*% qq
    }
    
    rr <- ww
    if (a>1)
      for (j in 1:(a-1))
        rr <- rr - (PP[,j] %*% ww) * RR[,j]
    
    tt <- scale(X %*% rr, scale=FALSE)   # centered X scores
    tnorm <- sqrt(sum(tt*tt))
    tt <- tt / tnorm
    rr <- rr / tnorm
    pp <- crossprod(X, tt)           # X loadings
    qq <- crossprod(rr, XtY)         # Y loadings, row vector!
    uu <- Y %*% matrix(qq, ncol=1)   # Y block factor scores
    if (a > 1)
      uu <- uu - TT %*% crossprod(TT, uu) # uu orth to previous tt values
    
    ## Now deflate crossprod matrices
    XtY <- XtY - (pp %*% qq)
    
    ## store weights and loadings
    TT[,a] <- tt
    PP[,a] <- pp 
    QQ[,a] <- qq
    RR[,a] <- rr
    UU[,a] <- uu

    if (!is.na(i <- match(a, ncomp))) {
      B[ , , i] <- RR[, 1:a, drop=FALSE] %*% t(QQ[, 1:a, drop=FALSE])
      if (!is.null(newX))
        Ypred[ , , i] <- newX %*% B[ , , i]
    }
  }

  if (!is.null(newX))
    list(B=B, Ypred=Ypred, Xscores=TT, Xload=PP, Yload=QQ, Yscores=UU)
  else
    list(B=B, Xload=PP, Yload=QQ, Xscores=TT, Yscores=UU)
}


