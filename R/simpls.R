### Contains R/Splus code for SIMPLS.
### Adapted from: De Jong, Chemolab v18 (1993) 251-263
### input: centred X and Y matrices, and the numbers of latent variables.

simpls <- function(X, Y, ncomp, newX=NULL)
{
  nobj <- dim(X)[1] # n in paper
  nvar <- dim(X)[2] # p in paper
  npred <- dim(Y)[2]

  S <- crossprod(X, Y)
  RR <- matrix(0, ncol=max(ncomp), nrow=nvar)
  PP <- matrix(0, ncol=max(ncomp), nrow=nvar)
  QQ <- matrix(0, ncol=max(ncomp), nrow=npred)
  TT <- matrix(0, ncol=max(ncomp), nrow=nobj)
  VV <- matrix(0, ncol=max(ncomp), nrow=nvar)
  B <- array(0, c(dim(X)[2], dim(Y)[2], length(ncomp)))

  if (!is.null(newX))
    Ypred <- array(0, c(dim(newX)[1], npred, length(ncomp)))

  for (a in 1:max(ncomp)) {
    qq <- svd(S)$v[,1]	        # Y block factor weights
    rr <- S %*% qq		# X block factor weights
    tt <- scale(X %*% rr, scale=FALSE)	 # center scores
    tnorm <- sqrt(sum(tt*tt))
    tt <- tt/tnorm		# normalize scores
    rr <- rr/tnorm		# adapt weights accordingly
    pp <- crossprod(X, tt)	# X block factor loadings
    qq <- crossprod(Y, tt)	# Y block factor loadings
    uu <- Y %*% qq		# Y block factor scores
    vv <- pp			# init orthogonal loadings
    if (a > 1){
      vv <- vv - VV %*% crossprod(VV, pp) # vv orth to previous loadings
      uu <- uu - TT %*% crossprod(TT, uu) # uu orth to previous tt values
    }
    vv <- vv/sqrt(sum(vv*vv))	# normalize orthogonal loadings
    S <- S - vv %*% crossprod(vv, S) 	# deflate S
    
    RR[,a] <- rr
    TT[,a] <- tt
    PP[,a] <- pp
    QQ[,a] <- qq
    VV[,a] <- vv

    if (!is.na(i <- match(a, ncomp))) {
      B[ , , i] <- RR[,1:a,drop=FALSE] %*% t(QQ[,1:a,drop=FALSE])
      if (!is.null(newX))
        Ypred[ , , i] <- newX %*% B[ , , i]
    }
  }

  XvarExpl <- diag(crossprod(PP)) / (sum(diag(var(X))) * (nobj - 1))
#  totalYvarExpl <- diag(crossprod(QQ)) / (sum(diag(var(Y))) * (nobj - 1))
  YvarExpl <- matrix(0, length(ncomp), npred)
  for (i in 1:max(ncomp))
    YvarExpl[i,] <- diag(cor(Y, X %*% B[ , , i]))^2
 # browser()
      
  if (!is.null(newX))
    list(B=B, XvarExpl=matrix(cumsum(XvarExpl)[ncomp], ncol=1),
         YvarExpl=YvarExpl, Ypred=Ypred)
  else
    list(B=B, XvarExpl=matrix(cumsum(XvarExpl)[ncomp], ncol=1),
         YvarExpl=YvarExpl)
}


