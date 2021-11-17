####FALCO####
library(parallel)
library(irlba)

#' Perform factor analysis in data with dependent samples.
#' 
#' Uses an algorithm called FALCO to estimate the average variance multipliers, latent loadings L, factors C, and eigenvalues of 1/pL(C'C)L'.
#' 
#' @param Y A \code{p} x \code{n} (number genes/methylation sites x sample size) matrix of observed expression/methylation
#' @param K The number of latent factors (an integer >0). If left unspecified, it is estimated using \code{CBCV_plus}.
#' @param Z A \code{n} x \code{r} (sample size x number of covariates) matrix of nuisance covariates. This will typically be the intercept (i.e. a vector of all ones). The default value is no covariates.
#' @param B An \code{n} x \code{n} positive semidefinite matrix or list of \code{n} x \code{n} positive semidefinite matrices that determine the correlation structure between samples. For example, these can be kinship matrices, partition matrices that group samples by blocks, diagonal matrices with 0's or 1's along the diagonal grouping samples that have the same residual variance, etc. The software includes the identity if it is not in the span of the user-provided matrices. The default is to treat all samples as independent with the same residual variance.
#' @param updateDelta \code{T} or \code{F}. If \code{T}, then FALCO performs one iteration PCA assuming Y ~ (LC', Delta, V), where Delta is a diagonal matrix. This tends to improve estimates for C if genomic units have very heterogeneous variances. However, given the computational cost, it is recommended this be set to \code{F} when p >> 10^4 (i.e. DNA methylation data). Defaults to  \code{F}.
#' @param A.ine A #inequality constaints x b matrix, where b = length(B) = num. variance components. \code{A.ine \%*\% tau >= 0}, where tau = (tau_1^2,...,tau_b^2) are the variance multipliers. The default is all variance multipliers are greater than 0.
#' @param A.equ A #equality constraints x b matrix, where \code{A.equ \%*\% tau = 0}. Defaults to no equality constraints.
#' @param c.ine A #inequality constraints vector, where \code{A.ine \%*\% tau - c.ine >= 0}. It is highly recommended users do NOT input anything other than \code{0}. Defaults to \code{0}.

#' @return A list \item{B}{The list of matrices used to parametrize the correlation between samples} \item{v}{The variance multipliers for each B, where v_1B_1 + ... + v_bB_b is the average variance (across genes) between samples} \item{Delta}{If updateDelta=\code{T}, the estimate for Delta in the model Y ~ (LC', Delta, V). If updateDelta=\code{F}, this is \code{NULL}.} \item{Lambda}{Estimates for the non-zero eigenvalues of 1/pL(C'C)L' if updateDelta=\code{F}. If updateDelta=\code{T}, returns estimates for the non-zero eigenvalues of the variance-standardized latent factors, 1/pDelta^\{-1/2\}L(C'C)L'Delta^\{-1/2\}.} \item{C}{The estimate for the latent factors. This has orthogonal columns, and is orthogonal to Z} \item{L}{The estimate for the latent loadings.}
#' @export
FALCO <- function(Y, K=NULL, Z=NULL, B=NULL, updateDelta=F, A.ine=NULL, A.equ=NULL, c.ine=NULL) {
  if (!is.null(B)) {
    if (is.matrix(B)) {
      B <- list(B)
    }
    B <- IncludeIdent(B)
  }
  if (is.null(K)) {
    cat("Using CBCV+ with nFolds=5 to estimate K...")
    tmp.K <- CBCV_plus(Y = Y, Cov = Z, maxK = ncol(Y)/2, B = B, nFolds = 5, simpleDelta = T, A.ine = A.ine, c.ine = c.ine, A.equ = A.equ, svd.method = "slow", plotit = T)
    K <- tmp.K$K.hat
    cat("done!\n")
  }
  if (is.null(B)) {
    cat(paste0("B is undefined. Either define B, or consider using standard PCA with ", K, " factors.\n", collapse = ""))
    return(list(K=K))
  }
  if (K < 0) {
    stop("The number of factors, K, is less than 0")
  }
  if (K == 0) {
    cat("The number of factors, K, is 0. No factor analysis needed.\n")
    return(0)
  }
  cat(paste0("FALCO with K=", K, " factors..."))
  
  B.orig <- B
  if (!is.null(Z)) {
    if (NROW(Z) != ncol(Y)) {stop("nrow(Z) must be the same as ncol(Y)")}
    Q.Z <- Compute.QQ(Z)
    Y <- Y%*%Q.Z
    B <- lapply(B, function(x){t(Q.Z)%*%x%*%Q.Z})
  }
  out <- list()
  out$B <- B.orig
  out$Delta <- NULL
  out.C <- EstimateCperp(Y=Y, K=K, X=NULL, Z=NULL, B=B, simpleDelta=T, A.ine=A.ine, c.ine=c.ine, A.equ=A.equ, svd.method="slow", return.all=F)
  v.hat <- out.C$rho
  
  p <- nrow(Y); n <- ncol(Y)
  V <- CreateV(B, v.hat)
  sqrt.V <- sqrt.mat2(V)
  svd.Ytilde <- svd(1/sqrt(p)*Y%*%sqrt.V$Rinv, nu = K, nv = K)
  C.hat <- cbind(sqrt.V$R%*%(svd.Ytilde$v[,1:K]))   #Need to multiply this by sqrt(n)
  tmp.v <- Est.delta(Y = Y, Cov=C.hat, V = sqrt.V, is.sqrt=T)
  
  if (updateDelta) {
    cat("\nRefining estimate for C and L...")
    out <- Refine.C(Y = Y, Z = NULL, C = C.hat, B = B, v = v.hat*tmp.v, A.ine = A.ine, rotate.C = T)
    if (!is.null(Z)) {out$C <- Q.Z%*%out$C}
    out$B <- B.orig
    cat("done!\n")
    return(out)
  }
  
  gamma.hat <- svd.Ytilde$d[1:K]^2
  L.hat <- cbind(svd.Ytilde$u[,1:K])%*%diag(sqrt(gamma.hat),nrow=K,ncol=K)    #Need to divide this by sqrt(n)
  R <- sqrt.mat2(t(C.hat)%*%C.hat)$R
  
  tmp.eigen <- eigen(R%*%diag(gamma.hat - tmp.v, nrow=K, ncol=K)%*%R, symmetric = T)
  U <- tmp.eigen$vectors
  out$v <- v.hat*tmp.v
  out$Lambda <- tmp.eigen$values
  out$C <- sqrt(n)*C.hat%*%solve(R)%*%U
  out$L <- 1/sqrt(n)*L.hat%*%R%*%U
  if (!is.null(Z)) {out$C <- Q.Z%*%out$C}
  cat("done!\n")
  return(out)
}

Refine.C <- function(Y, Z=NULL, C, B, v, A.ine=NULL, rotate.C=T) {
  Q <- Compute.QQ(cbind(Z,C))
  tmp <- Est.Corr.multB(Y = Y%*%Q, B = lapply(B,function(x){t(Q)%*%x%*%Q}), theta.0 = v, simple.rho = F, A = A.ine)
  out <- list(Delta=tmp$Delta, v=tmp$Rho)
  K <- ncol(cbind(C))
  p <- nrow(Y)
  
  if (!is.null(Z)) {
    Q.Z <- Compute.QQ(Z)
    n <- ncol(Q.Z)
    V <- CreateV(B = lapply(B,function(x){t(Q.Z)%*%x%*%Q.Z}), Rho = tmp$Rho)
    tmp.V <- sqrt.mat2(V)
    s <- svd( (Y%*%(Q.Z%*%tmp.V$Rinv))*sqrt(1/tmp$Delta), nu = K, nv = K )
    out$C <- Q.Z%*%tmp.V$R%*%cbind(s$v[,1:K])
  } else {
    n <- ncol(Y)
    V <- CreateV(B = B, Rho = tmp$Rho)
    tmp.V <- sqrt.mat2(V)
    s <- svd( (Y%*%tmp.V$Rinv)*sqrt(1/tmp$Delta), nu = K, nv = K )
    out$C <- tmp.V$R%*%cbind(s$v[,1:K])
  }

  out$L <- (cbind(s$u[,1:K])%*%diag(s$d[1:K])) * sqrt(tmp$Delta)
  if (!rotate.C || K==1) {
    if (K==1) {
      scale.C <- 1/sqrt(sum(out$C^2)/n)
      out$C <- out$C*scale.C
      out$L <- out$L/scale.C
    }
    return(out)
  }
  sqrt.CtC <- sqrt.mat2(1/n*t(out$C)%*%out$C)
  eig <- eigen(x = sqrt.CtC$R%*%diag(s$d[1:K]^2/p-1)%*%sqrt.CtC$R, symmetric = T)
  out$C <- out$C%*%(sqrt.CtC$Rinv%*%eig$vectors)
  out$L <- out$L%*%(sqrt.CtC$R%*%eig$vectors)
  out$Lambda <- eig$values*n
  return(out)
}

Compute.QQ <- function(X) {
  X <- cbind(X)
  qr.X <- qr(X)
  return( qr.Q(qr.X, complete = T)[,(qr.X$rank+1):nrow(X)] )
}

Est.delta <- function(Y, Cov=NULL, V, is.sqrt=F) {
  if (is.sqrt) {
    tmp <- 1/2*(V$Rinv + t(V$Rinv))
  } else {
    tmp <- sqrt.mat2(V)$Rinv; tmp <- 1/2*(tmp + t(tmp))
  }
  if (!is.null(Cov)) {
    Cov <- cbind(Cov)
    Y <- Y %*% tmp; Cov <- tmp %*% Cov
    R <- Y - Y %*% Cov %*% solve(t(Cov)%*%Cov,t(Cov))
    return(sum(R^2) / nrow(Y) / (nrow(Cov)-ncol(Cov)))
  }
  return( sum((Y %*% tmp)^2) / nrow(Y) / ncol(Y) )
}