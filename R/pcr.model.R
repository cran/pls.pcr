pcr.model <- function(X, Y, ncomp, newX)
{
  B <- array(0, c(dim(X)[2], dim(Y)[2], length(ncomp)))
  if (!is.null(newX))
    Ypred <- array(0, c(dim(newX)[1], dim(Y)[2], length(ncomp)))
  
  huhn <- La.svd(X)
  U <- huhn$u[,1:max(ncomp),drop=FALSE]
  D <- huhn$d[1:max(ncomp)]
  Vt <- huhn$vt[1:max(ncomp),,drop=FALSE]

  yve <- matrix(0, length(ncomp), dim(Y)[2])

  for (i in seq(along=ncomp)) {
    B[ , , i] <- t(Vt[1:ncomp[i], , drop=FALSE]) %*%
      diag(1/D[1:ncomp[i]], nrow=ncomp[i]) %*%
        t(U[ ,1:ncomp[i], drop=FALSE]) %*% Y
    if (!is.null(newX))
      Ypred[ , , i] <- newX %*% B[ , , i]

    yve[i,] <- diag(cor(Y, X %*% B[ , , i]))^2
  }

  if (!is.null(newX))
    list(B=B, XvarExpl=matrix(cumsum(D^2)[ncomp]/sum(huhn$d^2), ncol=1),
         YvarExpl=yve, Ypred=Ypred)
  else
    list(B=B, XvarExpl=matrix(cumsum(D^2)[ncomp]/sum(huhn$d^2), ncol=1),
         YvarExpl=yve)
}
