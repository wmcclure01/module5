#Wayne McClure Module #5 Doing Math LIS4370

##make the Matrices
A=matrix(1:100, nrow=10)
B=matrix(1:1000, nrow=10)

detA <- det(A)

#function to apply MP-pseudoinverse
exp.mat<-function(MAT, EXP, tol=NULL){
  MAT <- as.matrix(MAT)
  matdim <- dim(MAT)
  if(is.null(tol)){
    tol=min(1e-7, .Machine$double.eps*max(matdim)*max(MAT))
  }
  if(matdim[1]>=matdim[2]){ 
    svd1 <- svd(MAT)
    keep <- which(svd1$d > tol)
    res <- t(svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep]))
  }
  if(matdim[1]<matdim[2]){ 
    svd1 <- svd(t(MAT))
    keep <- which(svd1$d > tol)
    res <- svd1$u[,keep]%*%diag(svd1$d[keep]^EXP, nrow=length(keep))%*%t(svd1$v[,keep])
  }
  return(res)
}

invB <- exp.mat(B, -1)


head.matrix(invB)

library(matlib)
detA <- det(A)
detA

#invA <- inv(A)
svdA <- SVD(A)
svdA$v%*%solve(diag(svdA$d))%*%t(svdA$u)
solve(svdA$u %*% diag(svdA$d) %*% t(svdA$v))
