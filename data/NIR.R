"NIR" <- 
  list(matrix(scan("xtrainNIR",quiet=TRUE), nrow=21, byrow=TRUE),
       matrix(scan("xtestNIR",quiet=TRUE), nrow=7, byrow=TRUE),
       scan("ytrainNIR", quiet=TRUE),
       scan("ytestNIR", quiet=TRUE))
names(NIR) <- c("Xtrain","Xtest","Ytrain","Ytest")
