library(Rcpp)
library(MTS)
library(mvtnorm)
library(mc2d)
library(MASS)
library(L1pack)
library(CompQuadForm)
library(future.apply)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
Rcpp::sourceCpp("D:/R/code/self code/varmalikelihood.cpp")



check_and_fix_beta0 <- function(beta0, d, p, q) {#检查初始值 beta0 是否平稳可逆
  d2 <- d * d
  max_eval_PH <- 0
  max_eval_TH <- 0
  
  if (p > 0) {
    PH_mat <- matrix(beta0[1:(p * d2)], nrow = d, ncol = d * p)
    comp_PH <- matrix(0, d * p, d * p)
    comp_PH[1:d, ] <- PH_mat
    if (p > 1) {
      comp_PH[(d + 1):(d * p), 1:(d * (p - 1))] <- diag(d * (p - 1))
    }
    max_eval_PH <- max(Mod(eigen(comp_PH, only.values = TRUE)$values))
  }
  
  if (q > 0) {
    TH_mat <- matrix(beta0[(p * d2 + 1):((p + q) * d2)], nrow = d, ncol = d * q)
    comp_TH <- matrix(0, d * q, d * q)
    comp_TH[1:d, ] <- TH_mat
    if (q > 1) {
      comp_TH[(d + 1):(d * q), 1:(d * (q - 1))] <- diag(d * (q - 1))
    }
    max_eval_TH <- max(Mod(eigen(comp_TH, only.values = TRUE)$values))
  }
  
  if (max_eval_PH >= 0.999 || max_eval_TH >= 0.999) {
    warning("Initial beta0 violates stationarity/invertibility. Resetting to 0.")
    return(rep(0, length(beta0)))
  }
  return(beta0)
}


zzr.varma.est <- function(zt, beta0, p, q) {#VARMA(p,q) 估计
  d <- nrow(zt)
  beta0 <- check_and_fix_beta0(beta0, d, p, q)
  
  last_beta <- NULL      
  cached_result <- NULL  
  
  update_cache <- function(beta) {
    if (is.null(last_beta) || !identical(beta, last_beta)) {
      cached_result <<- zzr_varma_adjoint_list_rcpp(zt, beta, p, q)
      last_beta <<- beta
    }
  }
  
  target_fn <- function(beta) { update_cache(beta); return(cached_result$objective) }
  target_gr <- function(beta) { update_cache(beta); return(cached_result$gradient) }
  
  fit <- optim(
    par = beta0,
    fn = target_fn,
    gr = target_gr,
    method = "BFGS",
    control = list(maxit = 1000, trace = 0)
  )
  
  return(fit$par)
}


zzr.varma.boot.est <- function(zt, beta0, w, p, q) {#Bootstrap 加权 VARMA(p,q) 估计
  d <- nrow(zt)
  beta0 <- check_and_fix_beta0(beta0, d, p, q)
  
  last_beta <- NULL      
  cached_result <- NULL  
  
  update_cache <- function(beta) {
    if (is.null(last_beta) || !identical(beta, last_beta)) {
      cached_result <<- zzr_varma_boot_adjoint_list_rcpp(zt, beta, w, p, q)
      last_beta <<- beta
    }
  }
  
  target_fn <- function(beta) { update_cache(beta); return(cached_result$objective) }
  target_gr <- function(beta) { update_cache(beta); return(cached_result$gradient) }
  
  fit <- optim(
    par = beta0,
    fn = target_fn,
    gr = target_gr,
    method = "BFGS",
    control = list(maxit = 1000, trace = 0)
  )
  
  return(fit$par)
}



###########################################################################################
# 分块托普利兹矩阵构造
ToeplitzBlock <- function(X, lag.max) {
  np <- dim(X)
  n <- np[2]
  p <- np[1]
  X <- t(X)
  m <- lag.max + 1
  

  Accmat <- stats::acf(X, lag.max = lag.max, plot = FALSE, type = "covariance")$acf
  inveseC0 <- solve(Accmat[1,,])
  L <- t(chol(inveseC0))
  Lt <- t(L) 
  

  unique_blocks_lower <- vector("list", m)
  unique_blocks_upper <- vector("list", m)
  
  for (k in 1:m) {
    Ak <- Accmat[k,,]

    unique_blocks_lower[[k]] <- Lt %*% t(Ak) %*% L  
    unique_blocks_upper[[k]] <- Lt %*% Ak %*% L
  }
  

  out <- matrix(0, nrow = p * m, ncol = p * m)
  
  for (i in 0:lag.max) {
    for (j in i:lag.max) {
      lag_idx <- (j - i) + 1  
      
      row_idx <- (j * p + 1):(p * (j + 1))
      col_idx <- (i * p + 1):(p * (i + 1))
      

      out[row_idx, col_idx] <- unique_blocks_lower[[lag_idx]]
      if (i != j) {
        out[col_idx, row_idx] <- unique_blocks_upper[[lag_idx]]
      }
    }
  }
  return(out)
}

# 计算统计量
MahdiMcLeod <- function(X, lag.max) {
  np <- dim(X)
  n <- np[2]
  d <- np[1]
  

  mat <- ToeplitzBlock(X = X, lag.max = lag.max)
  

  det_res <- determinant(mat, logarithm = TRUE)
  Det1 <- (-n) * as.numeric(det_res$modulus)

  return(Det1)
}


#####################白噪声检验
MahdiMcLeod.wntest.boot1<-function(X,k_max,boot=1000){#白噪声检验, 权重为标准指数分布Exp(1)
  d = dim(X)[1]
  n = dim(X)[2]
  Y<-X
  X<-X-apply(X,1,mean)
  G0=solve(X%*%t(X)/n)

  ll=rep(0,k_max)
  ls=rep(0,k_max)
  rho.vec=rep(0,(d^2)*k_max)
  for (k in (1:k_max)){
    sm=matrix(0,d,d)
    X1<-X[,(k+1):n]
    X2<-X[,1:(n-k)]
    sm<-X1%*%t(X2)/n
    start<-(k-1)*(d^2)+1
    end<-k*(d^2)
    rho.vec[start:end]<-c(sm)
  }
  Q=MahdiMcLeod(X=X,lag.max=k_max)
  
  rho.vec.boot<-matrix(0,(d^2)*k_max,boot)
  sigma.vec.boot<-matrix(0,(d^2),boot)
  for(b in 1:boot){
    w<-rexp(n, rate = 1)
    Y.boot<-t(t(Y)*(matrix(rep(w, d), n, d)))
    X.boot<-Y-apply(Y.boot,1,mean)

    sigma.boot<-X.boot%*%t(X.boot)/n
    sigma.vec.boot[,b]<-c(sigma.boot)
    X.boot2<-t(t(X.boot)*(matrix(rep(w, d), n, d)))
    for(k in 1:k_max){
      X1.boot<-X.boot2[,(k+1):n]
      X2.boot<-X.boot[,1:(n-k)]
      sm.boot<-X1.boot%*%t(X2.boot)/n
      start<-(k-1)*(d^2)+1
      end<-k*(d^2)
      rho.vec.boot[start:end,b]<-c(sm.boot)
    }
  }
  cv<-sqrt(n)*(rho.vec.boot-rho.vec)
  V.boot<-cov(t(cv))
  sigma.mean<-matrix(apply(sigma.vec.boot,1,mean),d,d)
  svd.s<-svd(sigma.mean)
  sinv.root<-svd.s$u %*% diag(1/sqrt(svd.s$d)) %*% t(svd.s$v)
  K<-kronecker(diag(k_max),sinv.root)
  K<-kronecker(K,sinv.root)
  
  V<-K%*%V.boot%*%K
  
  
  
  M<-matrix(0,(k_max*d^2),(k_max*d^2))

  for(k in 1:k_max){
    start<-(k-1)*(d^2)+1
    end<-k*(d^2)
    M[start:end,start:end]<-sqrt(k_max-k+1)*diag(d^2)
    
  }
  V<-M%*%V%*%M
  lambda<-eigen(V)$values
  p.value<-davies(q=Q,lambda=lambda)$Qq#p值
  reject5<-(p.value<=0.05)###基于5%显著性水平是否拒绝假设
  reject10<-(p.value<=0.1)###基于10%显著性水平是否拒绝假设
  return(list(
    statistic = Q,           # 检验统计量
    p_value = p.value,          # p值
    reject_5pct = reject5,      # 5% 显著性水平拒绝结果 (TRUE/FALSE)
    reject_10pct = reject10     # 10% 显著性水平拒绝结果 (TRUE/FALSE)
  ))

}



MahdiMcLeod.wntest.boot2<-function(X,k_max,boot=1000){#白噪声检验, 权重为Mammen分布
  d = dim(X)[1]
  n = dim(X)[2]
  Y<-X
  X<-X-apply(X,1,mean)
  G0=solve(X%*%t(X)/n)
  
  ll=rep(0,k_max)
  ls=rep(0,k_max)
  rho.vec=rep(0,(d^2)*k_max)
  for (k in (1:k_max)){
    sm=matrix(0,d,d)
    X1<-X[,(k+1):n]
    X2<-X[,1:(n-k)]
    sm<-X1%*%t(X2)/n
    start<-(k-1)*(d^2)+1
    end<-k*(d^2)
    rho.vec[start:end]<-c(sm)
  }
  Q=MahdiMcLeod(X=X,lag.max=k_max)
  
  rho.vec.boot<-matrix(0,(d^2)*k_max,boot)
  sigma.vec.boot<-matrix(0,(d^2),boot)
  probuse<-(1+sqrt(5))/(2*sqrt(5))
  for(b in 1:boot){
    w<-rbern(n, prob=probuse)
    c<-which(w==1)
    w[c]=(3-sqrt(5))/2
    w[-c]=(3+sqrt(5))/2
    Y.boot<-t(t(Y)*(matrix(rep(w, d), n, d)))
    X.boot<-Y-apply(Y.boot,1,mean)
    
    sigma.boot<-X.boot%*%t(X.boot)/n
    sigma.vec.boot[,b]<-c(sigma.boot)
    X.boot2<-t(t(X.boot)*(matrix(rep(w, d), n, d)))
    for(k in 1:k_max){
      X1.boot<-X.boot2[,(k+1):n]
      X2.boot<-X.boot[,1:(n-k)]
      sm.boot<-X1.boot%*%t(X2.boot)/n
      start<-(k-1)*(d^2)+1
      end<-k*(d^2)
      rho.vec.boot[start:end,b]<-c(sm.boot)
    }
  }
  cv<-sqrt(n)*(rho.vec.boot-rho.vec)
  V.boot<-cov(t(cv))
  sigma.mean<-matrix(apply(sigma.vec.boot,1,mean),d,d)
  svd.s<-svd(sigma.mean)
  sinv.root<-svd.s$u %*% diag(1/sqrt(svd.s$d)) %*% t(svd.s$v)
  K<-kronecker(diag(k_max),sinv.root)
  K<-kronecker(K,sinv.root)
  
  V<-K%*%V.boot%*%K
  
  
  
  M<-matrix(0,(k_max*d^2),(k_max*d^2))
  
  for(k in 1:k_max){
    start<-(k-1)*(d^2)+1
    end<-k*(d^2)
    M[start:end,start:end]<-sqrt(k_max-k+1)*diag(d^2)
    
  }
  V<-M%*%V%*%M
  lambda<-eigen(V)$values
  p.value<-davies(q=Q,lambda=lambda)$Qq#p值
  reject5<-(p.value<=0.05)###基于5%显著性水平是否拒绝假设
  reject10<-(p.value<=0.1)###基于10%显著性水平是否拒绝假设
  
  return(list(
    statistic = Q,           # 检验统计量
    p_value = p.value,          # p值
    reject_5pct = reject5,      # 5% 显著性水平拒绝结果 (TRUE/FALSE)
    reject_10pct = reject10     # 10% 显著性水平拒绝结果 (TRUE/FALSE)
  ))
  
}


########################################################
zzr.varma.gof.boot1<-function(zt,beta0,p,q,k_max,boot=1000){#VARMA(p,q)拟合优度检验, 权重为标准指数分布Exp(1)
  nd<-dim(zt)
  n<-nd[2]
  d<-nd[1]
  
  theta<-zzr.varma.est(zt=zt,beta0=beta0,p=p,q=q)
  if(p==0){
    PH1<-matrix(0,d,d)
  }else{
    PH1<-matrix(theta[1:((d*d)*p)],d,(d*p))
  }
  
  if(q==0){
    TH1<-matrix(0,d,d)
  }else{
    TH1<-matrix(theta[((d*d)*p+1):(d*d*(p+q))],d,(d*q))
  }

  res<-zzrvarmaResiduals_cpp(zt=zt,PH=PH1,TH=TH1,p=p,q=q)
  Y<-res
  X<-res
  G0=solve(X%*%t(X)/n)

  ll=rep(0,k_max)
  ls=rep(0,k_max)
  rho.vec <- rep(0, (d^2) * k_max)
  

  for (k in 1:k_max) {
    X1 <- X[, (k+1):n, drop=FALSE]
    X2 <- X[, 1:(n-k), drop=FALSE]
    sm <- tcrossprod(X1, X2) / n   # 等价于 X1 %*% t(X2) / n
    
    start <- (k-1)*(d^2) + 1
    end <- k*(d^2)
    rho.vec[start:end] <- c(sm)
  }
  
  Q=MahdiMcLeod(X=X,lag.max=k_max)

  
  
  rho.vec.boot<-matrix(0,(d^2)*k_max,boot)
  sigma.vec.boot<-matrix(0,(d^2),boot)
  #probuse<-(1+sqrt(5))/(2*sqrt(5))
  for(b in 1:boot){
    w<-rexp(n=n,rate=1)
    #w<-rbern(n, prob=p)
    #c<-which(w==1)
    #w[c]=(3-sqrt(5))/2
    #w[-c]=(3+sqrt(5))/2
    theta.boot<-zzr.varma.boot.est(zt=zt,beta0=beta0,w=w,p=p,q=q)
    
    if(p==0){
      PH.boot<-matrix(0,d,d)
    }else{
      PH.boot<-matrix(theta.boot[1:((d*d)*p)],d,(d*p))
    }
    
    if(q==0){
      TH.boot<-matrix(0,d,d)
    }else{
      TH.boot<-matrix(theta.boot[((d*d)*p+1):(d*d*(p+q))],d,(d*q))
    }
    
    

    res.boot<-zzrvarmaResiduals_cpp(zt=zt,PH=PH.boot,TH=TH.boot,p=p,q=q)

    X.boot<-res.boot

    cpp_res <- zzr_boot_cov_cpp(X.boot, w, k_max)
    sigma.vec.boot[, b] <- cpp_res[1:(d^2)]
    rho.vec.boot[, b] <- cpp_res[(d^2 + 1):length(cpp_res)]
  }
  cv<-sqrt(n)*(rho.vec.boot-rho.vec)

  cv_mean <- rowMeans(cv)
  

  cv_centered <- cv - cv_mean
  

  V.boot <- tcrossprod(cv_centered) / boot 
  
  
  sigma.mean <- matrix(rowMeans(sigma.vec.boot), d, d)
  svd.s<-svd(sigma.mean)
  sinv.root<-svd.s$u %*% diag(1/sqrt(svd.s$d)) %*% t(svd.s$v)
  K<-kronecker(diag(k_max),sinv.root)
  K<-kronecker(K,sinv.root)
  
  V2<-K%*%V.boot%*%K

  m_diag <- rep(sqrt(k_max - 1:k_max + 1), each = d^2)
  V1 <- V2 * tcrossprod(m_diag)

  
  lambda1 <- eigen(V1, symmetric = TRUE, only.values = TRUE)$values
  p.value<-davies(q=Q,lambda=lambda1)$Qq
  reject5<-(p.value<=0.05)###基于5%显著性水平是否拒绝假设
  reject10<-(p.value<=0.1)###基于10%显著性水平是否拒绝假设
  
  return(list(
    statistic = Q,           # 检验统计量
    p_value = p.value,          # p值
    reject_5pct = reject5,      # 5% 显著性水平拒绝结果 (TRUE/FALSE)
    reject_10pct = reject10     # 10% 显著性水平拒绝结果 (TRUE/FALSE)
  ))
}




zzr.varma.gof.boot2<-function(zt,beta0,p,q,k_max,boot=1000){#VARMA(p,q)拟合优度检验, 权重为Mammen分布
  nd<-dim(zt)
  n<-nd[2]
  d<-nd[1]
  
  theta<-zzr.varma.est(zt=zt,beta0=beta0,p=p,q=q)
  if(p==0){
    PH1<-matrix(0,d,d)
  }else{
    PH1<-matrix(theta[1:((d*d)*p)],d,(d*p))
  }
  
  if(q==0){
    TH1<-matrix(0,d,d)
  }else{
    TH1<-matrix(theta[((d*d)*p+1):(d*d*(p+q))],d,(d*q))
  }
  
  res<-zzrvarmaResiduals_cpp(zt=zt,PH=PH1,TH=TH1,p=p,q=q)
  Y<-res
  X<-res
  G0=solve(X%*%t(X)/n)
  
  ll=rep(0,k_max)
  ls=rep(0,k_max)
  rho.vec <- rep(0, (d^2) * k_max)
  
  
  for (k in 1:k_max) {
    X1 <- X[, (k+1):n, drop=FALSE]
    X2 <- X[, 1:(n-k), drop=FALSE]
    sm <- tcrossprod(X1, X2) / n   # 等价于 X1 %*% t(X2) / n
    
    start <- (k-1)*(d^2) + 1
    end <- k*(d^2)
    rho.vec[start:end] <- c(sm)
  }
  
  Q=MahdiMcLeod(X=X,lag.max=k_max)
  
  
  
  rho.vec.boot<-matrix(0,(d^2)*k_max,boot)
  sigma.vec.boot<-matrix(0,(d^2),boot)
  probuse<-(1+sqrt(5))/(2*sqrt(5))
  for(b in 1:boot){
    #w<-rexp(n=n,rate=1)
    w<-rbern(n, prob=probuse)
    c<-which(w==1)
    w[c]=(3-sqrt(5))/2
    w[-c]=(3+sqrt(5))/2
    theta.boot<-zzr.varma.boot.est(zt=zt,beta0=beta0,w=w,p=p,q=q)
    
    if(p==0){
      PH.boot<-matrix(0,d,d)
    }else{
      PH.boot<-matrix(theta.boot[1:((d*d)*p)],d,(d*p))
    }
    
    if(q==0){
      TH.boot<-matrix(0,d,d)
    }else{
      TH.boot<-matrix(theta.boot[((d*d)*p+1):(d*d*(p+q))],d,(d*q))
    }
    
    
    
    res.boot<-zzrvarmaResiduals_cpp(zt=zt,PH=PH.boot,TH=TH.boot,p=p,q=q)
    
    X.boot<-res.boot
    
    cpp_res <- zzr_boot_cov_cpp(X.boot, w, k_max)
    sigma.vec.boot[, b] <- cpp_res[1:(d^2)]
    rho.vec.boot[, b] <- cpp_res[(d^2 + 1):length(cpp_res)]
  }
  cv<-sqrt(n)*(rho.vec.boot-rho.vec)
  
  cv_mean <- rowMeans(cv)
  
  
  cv_centered <- cv - cv_mean
  
  
  V.boot <- tcrossprod(cv_centered) / boot 
  
  
  sigma.mean <- matrix(rowMeans(sigma.vec.boot), d, d)
  svd.s<-svd(sigma.mean)
  sinv.root<-svd.s$u %*% diag(1/sqrt(svd.s$d)) %*% t(svd.s$v)
  K<-kronecker(diag(k_max),sinv.root)
  K<-kronecker(K,sinv.root)
  
  V2<-K%*%V.boot%*%K
  
  m_diag <- rep(sqrt(k_max - 1:k_max + 1), each = d^2)
  V1 <- V2 * tcrossprod(m_diag)
  
  
  lambda1 <- eigen(V1, symmetric = TRUE, only.values = TRUE)$values
  p.value<-davies(q=Q,lambda=lambda1)$Qq
  reject5<-(p.value<=0.05)###基于5%显著性水平是否拒绝假设
  reject10<-(p.value<=0.1)###基于10%显著性水平是否拒绝假设
  
  return(list(
    statistic = Q,           # 检验统计量
    p_value = p.value,          # p值
    reject_5pct = reject5,      # 5% 显著性水平拒绝结果 (TRUE/FALSE)
    reject_10pct = reject10     # 10% 显著性水平拒绝结果 (TRUE/FALSE)
  ))
}












##############



model1<-function(n,d,S){
  sigma=S
  Z<-t(mvrnorm((n+100),mu=rep(0,d),Sigma=sigma))
  X<-Z
  Y<-X[,101:(n+100)]
  #  Y<-Z[,101:(n+100)]
  return(Y)
}



sim1<-function(n,d,k_max,S){
  X<-model1(n=n,d=d,S=S)

  t1<-MahdiMcLeod.wntest.boot1(X=X,k_max=k_max,boot=1000)
  t2<-MahdiMcLeod.wntest.boot2(X=X,k_max=k_max,boot=1000)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))

}


d=4

r<-ceiling(d/2.5)
S<-matrix(0,d,d)
for(q in 1:(floor(d/r))){
  S[(r*(q-1)+1):(r*q),(r*(q-1)+1):(r*q)]=0.8
}
for(k in 1:d){
  S[k,k]=1
}




t1<- replicate(1000,sim1(n=500,d=d,k_max=5,S=S))
apply(t1,1,mean)

#######################################
model2<-function(n,d,C){
  X<-matrix(0,(n+30),d)
  S0<-matrix(0,d,d)
  S1<-matrix(0,d,d)
  for(i in 2:(n+30)){

    S1<-diag(d)+t(C)%*%X[i-1,]%*%t(X[i-1,])%*%C
    X[i,]<-MASS::mvrnorm(1,rep(0,d),S1)

  }
  Y<-t(X[31:(n+30),])
  return(Y)
}

sim2<-function(n,d,k_max,C){
  X<-model2(n=n,d=d,C=C)
  
  t1<-MahdiMcLeod.wntest.boot1(X=X,k_max=k_max,boot=1000)
  t2<-MahdiMcLeod.wntest.boot2(X=X,k_max=k_max,boot=1000)
  
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


d=4
vec.C<-rep(0,(d^2))
C1<-sort(sample((1:(d^2)),(2*d)))
vec.C[C1]<-runif((2*d),0.1,0.9)
C<-matrix(vec.C,nrow=d,ncol=d)


t2<- replicate(1000,sim2(n=500,d=d,k_max=5,C=C))
apply(t2,1,mean)

############
model3<-function(n,phi,theta,sig0){
  e<-t(mvrnorm((n+30),rep(0,2),sig0))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}


sim3<-function(n,k_max,PH,TH,sig0,beta0){

  zt<-model3(n=n,phi=PH,theta=TH,sig0=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0=matrix(c(1,0.9,0.9,1),2,2)

th1=matrix(c(0.5,-0.4,0,0.4),2,2)

beta0<-c(p1,th1)


t3<- replicate(1000,sim3(n=500,k_max=5,PH=p1,TH=th1,sig0=sig0,beta0=beta0))
apply(t3,1,mean)

###############
generate.garch<-function(n){
  error<-t(mvrnorm((n+30),rep(0,2),diag(2)))
  y<-matrix(0,2,(n+30))
  h1.old<-0
  h2.old<-0
  for(i in 2:(n+30)){
    h1.new<-0.6+0.3*y[1,i-1]^2+0.15*y[2,i-1]^2+0.25*h1.old
    h2.new<-0.4+0.5*y[2,i-1]^2+0.1*h1.old+0.1*h2.old
    
    y[1,i]<-sqrt(h1.new)*error[1,i]
    y[2,i]<-sqrt(h2.new)*error[2,i]
    
    h1.old<-h1.new
    h2.old<-h2.new
  }
  return(y[,31:(n+30)])
}


model4<-function(n,phi,theta){
  e<-generate.garch((n+30))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}


sim4<-function(n,k_max,PH,TH,beta0){

  zt<-model4(n=n,phi=PH,theta=TH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)

th1=matrix(c(0.5,-0.4,0,0.4),2,2)

beta0<-c(p1,th1)


t4<- replicate(1000,sim4(n=500,k_max=5,PH=p1,TH=th1,beta0=beta0))
apply(t4,1,mean)


##################
model5<-function(n,phi,theta,sigma){
  e<-t(rmLaplace((n+30),rep(0,2),sigma))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim5<-function(n,k_max,PH,TH,sig0,beta0){

  zt<-model5(n=n,phi=PH,theta=TH,sigma=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0=matrix(c(1,0.9,0.9,1),2,2)
th1=matrix(c(0.5,-0.4,0,0.4),2,2)

beta0<-c(p1,th1)

plan(multisession, workers = 12)
t5<- future_replicate(1000,sim5(n=500,k_max=5,PH=p1,TH=th1,sig0=sig0,beta0=beta0))
apply(t5,1,mean)

###############
generate.mix<-function(n){
  e<-t(rmvt((n+30),diag(2),df=5))
  error<-matrix(0,2,(n+30))
  for(i in 3:(n+30)){
    error[,i]<-e[,i]*e[,i-1]*e[,i-2]
  }
  return(error[,31:(n+30)])
}

model6<-function(n,phi,theta){
  e<-generate.mix((n+30))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim6<-function(n,k_max,PH,TH,beta0){

  zt<-model6(n=n,phi=PH,theta=TH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)

th1=matrix(c(0.5,-0.4,0,0.4),2,2)

beta0<-c(p1,th1)


t6<- replicate(1000,sim6(n=500,k_max=5,PH=p1,TH=th1,beta0=beta0))
apply(t6,1,mean)


######################################################################
model7<-function(n,phi,theta,sig0){
  e<-t(mvrnorm((n+30),rep(0,2),sig0))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]

  }
  return(X[,31:(n+30)])
}


sim7<-function(n,k_max,PH,TH,sig0,beta0){

  zt<-model7(n=n,phi=PH,theta=TH,sig0=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0=matrix(c(1,0.8,0.8,1),2,2)

th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1)


t7<- replicate(1000,sim7(n=500,k_max=5,PH=p1,TH=th1,sig0=sig0,beta0=beta0))
apply(t7,1,mean)

########################################################################
generate.garch<-function(n){
  error<-t(mvrnorm((n+30),rep(0,2),diag(2)))
  y<-matrix(0,2,(n+30))
  h1.old<-0
  h2.old<-0
  for(i in 2:(n+30)){
    h1.new<-0.6+0.3*y[1,i-1]^2+0.15*y[2,i-1]^2+0.25*h1.old
    h2.new<-0.4+0.5*y[2,i-1]^2+0.1*h1.old+0.1*h2.old
    
    y[1,i]<-sqrt(h1.new)*error[1,i]
    y[2,i]<-sqrt(h2.new)*error[2,i]
    
    h1.old<-h1.new
    h2.old<-h2.new
  }
  return(y[,31:(n+30)])
}


model8<-function(n,phi,theta){
  e<-generate.garch((n+30))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim8<-function(n,k_max,PH,TH,beta0){
  zt<-model8(n=n,phi=PH,theta=TH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)

th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1)


t8<- replicate(1000,sim8(n=500,k_max=5,PH=p1,TH=th1,beta0=beta0))
apply(t8,1,mean)


##############################
model9<-function(n,phi,theta,sigma){
  e<-t(rmLaplace((n+30),rep(0,2),sigma))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim9<-function(n,k_max,PH,TH,sig0,beta0){

  zt<-model9(n=n,phi=PH,theta=TH,sigma=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0=matrix(c(1,0.8,0.8,1),2,2)
th1=matrix(c(0.4,-0.3,0,0.2),2,2)
beta0<-c(p1)


t9<- replicate(1000,sim9(n=500,k_max=5,PH=p1,TH=th1,sig0=sig0,beta0=beta0))
apply(t9,1,mean)



################################
generate.mix<-function(n){
  e<-t(rmvt((n+30),diag(2),df=5))
  error<-matrix(0,2,(n+30))
  for(i in 3:(n+30)){
    error[,i]<-e[,i]*e[,i-1]*e[,i-2]
  }
  return(error[,31:(n+30)])
}

model10<-function(n,phi,theta){
  e<-generate.mix((n+30))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim10<-function(n,k_max,PH,TH,beta0){
  zt<-model10(n=n,phi=PH,theta=TH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)

th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1)


t10<- replicate(1000,sim10(n=500,k_max=5,PH=p1,TH=th1,beta0=beta0))
apply(t10,1,mean)




##############################
model12<-function(n,phi,sigma){
  e<-t(rmvt((n+30),sigma,df=4))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]
  }
  return(X[,31:(n+30)])
}

sim12<-function(n,k_max,PH,sig0,beta0){
  zt<-model12(n=n,phi=PH,sigma=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0<-matrix(c(1,0.8,0.8,1),2,2)


beta0<-c(p1)


t12<- replicate(1000,sim12(n=500,k_max=5,PH=p1,sig0=sig0,beta0=beta0))
apply(t12,1,mean)


###############################
model13<-function(n,phi,sigma){
  e<-t(rmvt((n+30),sigma,df=2))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]
  }
  return(X[,31:(n+30)])
}

sim13<-function(n,k_max,PH,sig0,beta0){

  zt<-model13(n=n,phi=PH,sigma=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0<-matrix(c(1,0.8,0.8,1),2,2)


beta0<-c(p1)


t13<- replicate(1000,sim13(n=500,k_max=5,PH=p1,sig0=sig0,beta0=beta0))
apply(t13,1,mean)

################################
generate.garch2<-function(n){
  error<-t(mvrnorm((n+30),rep(0,2),diag(2)))
  y<-matrix(0,2,(n+30))
  h1.old<-0
  h2.old<-0
  for(i in 2:(n+30)){
    h1.new<-0.6+0.4*y[1,i-1]^2+0.59*h1.old
    h2.new<-0.6+0.6*y[2,i-1]^2+0.39*h2.old
    
    y[1,i]<-sqrt(h1.new)*error[1,i]
    y[2,i]<-sqrt(h2.new)*error[2,i]
    
    h1.old<-h1.new
    h2.old<-h2.new
  }
  return(y[,31:(n+30)])
}


model14<-function(n,phi){
  e<-generate.garch2((n+30))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]
  }
  return(X[,31:(n+30)])
}


sim14<-function(n,k_max,PH,beta0){

  zt<-model14(n=n,phi=PH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)


beta0<-c(p1)


t14<- replicate(1000,sim14(n=500,k_max=5,PH=p1,beta0=beta0))
apply(t14,1,mean)

###############################
model15<-function(n,phi,theta,sigma){
  e<-t(rmvt((n+30),sigma,df=4))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim15<-function(n,k_max,PH,TH,sig0,beta0){

  zt<-model15(n=n,phi=PH,theta=TH,sigma=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0<-matrix(c(1,0.8,0.8,1),2,2)
th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1)

t15<- replicate(1000,sim15(n=500,k_max=5,PH=p1,TH=th1,sig0=sig0,beta0=beta0))
apply(t15,1,mean)


##############################
model16<-function(n,phi,theta,sigma){
  e<-t(rmvt((n+30),sigma,df=2))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}

sim16<-function(n,k_max,PH,TH,sig0,beta0){

  zt<-model16(n=n,phi=PH,theta=TH,sigma=sig0)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)
sig0<-matrix(c(1,0.8,0.8,1),2,2)
th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1)


t16<-replicate(1000,sim16(n=500,k_max=5,PH=p1,TH=th1,sig0=sig0,beta0=beta0))
apply(t16,1,mean)


######################################
generate.garch2<-function(n){
  error<-t(mvrnorm((n+30),rep(0,2),diag(2)))
  y<-matrix(0,2,(n+30))
  h1.old<-0
  h2.old<-0
  for(i in 2:(n+30)){
    h1.new<-0.6+0.4*y[1,i-1]^2+0.59*h1.old
    h2.new<-0.6+0.6*y[2,i-1]^2+0.39*h2.old
    
    y[1,i]<-sqrt(h1.new)*error[1,i]
    y[2,i]<-sqrt(h2.new)*error[2,i]
    
    h1.old<-h1.new
    h2.old<-h2.new
  }
  return(y[,31:(n+30)])
}


model17<-function(n,phi,theta){
  e<-generate.garch2((n+30))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}


sim17<-function(n,k_max,PH,TH,beta0){
  zt<-model17(n=n,phi=PH,theta=TH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=0,k_max=k_max)
  return(c(t1$reject_5pct,t1$reject_10pct,t2$reject_5pct,t2$reject_10pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)

th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1)


t17<- replicate(1000,sim17(n=500,k_max=5,PH=p1,TH=th1,beta0=beta0))
apply(t17,1,mean)


######################################
model18<-function(n,phi,theta){
  e<-t(mvrnorm((n+30),rep(0,2),diag(2)))
  X<-matrix(0,2,(n+30))
  for(i in 2:(n+30)){
    X[,i]<-phi%*%X[,i-1]+e[,i]-theta%*%e[,i-1]
  }
  return(X[,31:(n+30)])
}


sim18<-function(n,k_max,PH,TH,beta0){

  zt<-model18(n=n,phi=PH,theta=TH)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  return(c(t1$reject_5pct,t2$reject_5pct))
}

p1=matrix(c(0.2,-0.5,0.3,0),2,2)

th1=matrix(c(0.4,-0.3,0,0.2),2,2)

beta0<-c(p1,th1)


t18<- replicate(1000,sim18(n=500,k_max=5,PH=p1,TH=th1,beta0=beta0))
apply(t18,1,mean)



#############################
model19<-function(n,phi1,phi2,theta1,theta2){
  e<-t(mvrnorm((n+30),rep(0,2),diag(2)))
  X<-matrix(0,2,(n+30))
  for(i in 4:(n+30)){
    X[,i]<-phi1%*%X[,i-1]+phi2%*%X[,i-3]+e[,i]-theta1%*%e[,i-1]-theta2%*%e[,i-3]
  }
  return(X[,31:(n+30)])
}

sim19<-function(n,k_max,PH1,PH2,TH1,TH2,beta0){
  zt<-model19(n=n,phi1=PH1,phi2=PH2,theta1=TH1,theta2=TH2)
  t1<-zzr.varma.gof.boot1(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  t2<-zzr.varma.gof.boot2(zt=zt,beta0=beta0,p=1,q=1,k_max=k_max)
  return(c(t1$reject_5pct,t2$reject_5pct))
}


p1=matrix(c(0.2,-0.5,0.3,0),2,2)
p2=matrix(c(0,0.1,-0.1,0),2,2)
th1=matrix(c(0.4,-0.3,0,0.2),2,2)
th2=matrix(c(0.1,-0.1,0,0.05),2,2)
beta0<-c(p1,th1)



t19<- replicate(1000,sim19(n=500,k_max=5,PH1=p1,PH2=p2,TH1=th1,TH2=th2,beta0=beta0))
apply(t19,1,mean)


