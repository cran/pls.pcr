pcr.model <- function(X, Y, ncomp, newX)
{
  B <- array(0, c(dim(X)[2], dim(Y)[2], length(ncomp)))
  if (!is.null(newX))
    Ypred <- array(0, c(dim(newX)[1], dim(Y)[2], length(ncomp)))
  
  huhn <- La.svd(X)
  U <- huhn$u[,1:max(ncomp),drop=FALSE]
  D <- huhn$d[1:max(ncomp)]
  Vt <- huhn$vt[1:max(ncomp),,drop=FALSE]

  for (i in seq(along=ncomp)) {
    B[ , , i] <- t(Vt[1:ncomp[i], , drop=FALSE]) %*%
      diag(1/D[1:ncomp[i]], nrow=ncomp[i]) %*%
        t(U[ ,1:ncomp[i], drop=FALSE]) %*% Y
    if (!is.null(newX))
      Ypred[ , , i] <- newX %*% B[ , , i]
  }

  if (!is.null(newX))
    list(B=B, Ypred=Ypred,
         Xscores=U[,1:max(ncomp),drop=FALSE],
         Xload=t(diag(D) %*% Vt)[,1:max(ncomp),drop=FALSE])
  else
    list(B=B, Xscores=U[,1:max(ncomp),drop=FALSE],
         Xload=t(diag(D) %*% Vt)[,1:max(ncomp),drop=FALSE])
}
