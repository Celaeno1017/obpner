###estimate OBP under NER model############
library(lme4)
####input X: 3-D array with i: area i; j: sample j in n_i; k: covariate k in 1-p.
####input Y: 2-D array with i: arear i; j; sample j in n_i.
####input r: sampling fraction for each area i


#####fh area-level model##############
OBP_fh <- function (X, X_pop,Y, r,sigsq_e, sigsq_v){
  
X_bar <- X_pop
#X_bar <- apply(X,3,rowMeans)
Y_bar <- rowMeans(Y)  
n <- apply(X, 1, length)/dim(X_bar)[2]

Qr <- function(x,X_bar,Y_bar,n,r,sigsq_e){
  sigsq_v <- x[1]
  #sigsq_v <- x[2]
  
  tau <- diag((sigsq_e/n)/(sigsq_e/n+sigsq_v))
  D <- diag(sigsq_e/n)
  V <- tau^2
  xinv <- solve(t(X_bar)%*%V%*%X_bar)
  P <- diag(length(Y_bar))-tau%*%X_bar%*%xinv%*%t(X_bar)%*%tau
  return(sum(diag(D-2*tau%*%D))+(t(Y_bar)%*%tau%*%P%*%tau%*%Y_bar))
  
}

result <- optim(par = c(sigsq_v), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,n=n,r=r, sigsq_e=sigsq_e, lower = c(10^-9),method ="L-BFGS-B")

sigsq_v <- result$par

tau <- diag((sigsq_e/n)/(sigsq_e/n+sigsq_v))
beta <- solve(t(X_bar)%*%tau^2%*%X_bar)%*%t(X_bar)%*%tau^2%*%Y_bar
theta <- X_pop%*%beta + diag((sigsq_v)/(sigsq_e/n+sigsq_v))%*%(Y_bar-(X_pop)%*%beta)


return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))

return(result)
}


######fh unit-context model#############
OBP_fh_uc <- function (X, X_pop,Y, r,sigsq_e, sigsq_v){
  
  X_bar <- X_pop
  n <- apply(X, 1, length)/dim(X_bar)[2]
  n <- rep(1,sum(n))
  
  Y_unit <- as.vector(t(Y))
  X_unit <- kronecker(X_bar,rep(1,4))
  Qr <- function(x,X_bar,Y_bar,n,r,sigsq_e){
    sigsq_v <- x[1]
    #sigsq_v <- x[2]
    
    tau <- diag((sigsq_e/n)/(sigsq_e/n+sigsq_v))
    D <- diag(sigsq_e/n)
    V <- tau^2
    xinv <- solve(t(X_bar)%*%V%*%X_bar)
    P <- diag(length(Y_bar))-tau%*%X_bar%*%xinv%*%t(X_bar)%*%tau
    return(sum(diag(D-2*tau%*%D))+(t(Y_bar)%*%tau%*%P%*%tau%*%Y_bar))
    
  }
  
  result <- optim(par = c(sigsq_v), fn = Qr, X_bar = X_unit,Y_bar=Y_unit,n=n,r=r, sigsq_e=sigsq_e, lower = c(10^-9),method ="L-BFGS-B")
  
  sigsq_v <- result$par
  
  tau <- diag((sigsq_e/n)/(sigsq_e/n+sigsq_v))
  beta <- solve(t(X_unit)%*%tau^2%*%X_unit)%*%t(X_unit)%*%tau^2%*%Y_unit
  theta <- X_pop%*%beta + .colMeans(diag((sigsq_v)/(sigsq_e/n+sigsq_v))%*%(Y_unit-(X_unit)%*%beta),length(n)/length(apply(data$X_samp, 1, length)), length(apply(data$X_samp, 1, length)))
  
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
  
  return(result)
}


OBP_fh_model <- function (X, X_pop,Y, r,sigsq_e, sigsq_v){
  
  X_bar <- X_pop
  #X_bar <- apply(X,3,rowMeans)
  Y_bar <- rowMeans(Y)  
  n <- apply(X, 1, length)/dim(X_bar)[2]
  
  Qr <- function(x,X_bar,Y_bar,n,r,sigsq_e){
    sigsq_v <- x[1]
    #sigsq_v <- x[2]
    
    tau <- diag((1-r)*(sigsq_e)/(sigsq_e+n*sigsq_v))
    D <- diag(sigsq_e/n)
    V <- tau^2
    xinv <- solve(t(X_bar)%*%V%*%X_bar)
    P <- diag(length(Y_bar))-tau%*%X_bar%*%xinv%*%t(X_bar)%*%tau
    return(sum(diag(D-2*tau%*%D))+(t(Y_bar)%*%tau%*%P%*%tau%*%Y_bar))
    
  }
  
  result <- optim(par = c(sigsq_v), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,n=n,r=r, sigsq_e=sigsq_e, lower = c(10^-9),method ="L-BFGS-B")
  
  sigsq_v <- result$par
  
  # tau <- diag((1-r)*(sigsq_e)/(sigsq_e+n*sigsq_v))
  # beta <- solve(t(X_bar)%*%tau^2%*%X_bar)%*%t(X_bar)%*%tau^2%*%Y_bar
  # theta <- X_pop%*%beta + diag(r+(1-r)*(n*sigsq_v)/(sigsq_e+n*sigsq_v))%*%(Y_bar-(X_bar)%*%beta)
  
  X_bar <- apply(X,3,rowMeans)
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  H <- G%*%X_bar-X_pop
  Hinv <- solve(t(H)%*%H)
  
  
  beta <- Hinv%*%t(H)%*%(G-diag(length(Y_bar)))%*%Y_bar
  
  theta <- X_pop%*%beta+G%*%(Y_bar-X_bar%*%beta)
  
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
  
  return(result)
}



OBP_est_new <- function (X, X_pop, Y, r,sigsq_e, sigsq_v){
  X_bar <- apply(X,3,rowMeans)
  Y_bar <- rowMeans(Y)  
  n <- apply(X, 1, length)/dim(X_bar)[2]
  
  Qr <- function(x,X_bar,Y_bar,X_pop,n,r,sigsq_e){
    sigsq_v <- x[1]
    #sigsq_v <- x[2]
    
    tau <- diag((1-r)*(sigsq_e)/(sigsq_e+n*sigsq_v))
    D <- diag(sigsq_e/n)
    
    G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
    H <- G%*%X_bar-X_pop
    Hinv <- solve(t(H)%*%H)
    
    P <- diag(length(Y_bar))-H%*%Hinv%*%t(H)
    return(sum(diag(D-2*tau%*%D))+(t(Y_bar)%*%P%*%(G-diag(length(Y_bar)))^2%*%Y_bar))
    
  }
  
  
  result <- optim(par = c(sigsq_v), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,n=n,r=r, sigsq_e=sigsq_e, lower = c(sigsq_v*0.01+0.1),upper = c(sigsq_v*3000+0.2),method ="L-BFGS-B")
  
  sigsq_v <- result$par
  
  tau <- diag((1-r)*(sigsq_e)/(sigsq_e+n*sigsq_v))
  D <- diag(sigsq_e/n)
  
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  H <- G%*%X_bar-X_pop
  Hinv <- solve(t(H)%*%H)
  
  
  beta <- Hinv%*%t(H)%*%(G-diag(length(Y_bar)))%*%Y_bar
  theta <- X_pop%*%beta+G%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,beta=beta,theta=theta))
}

OBP_model <- function (X, X_pop, Y, r,sigsq_e, sigsq_v,beta){
  X_bar <- apply(X,3,rowMeans)
  Y_bar <- rowMeans(Y)  
  n <- apply(X, 1, length)/dim(X_bar)[2]
  
  Qr <- function(x,X_bar,Y_bar,X_pop,n,r){
    sigsq_v <- x[1]
    sigsq_e <- x[2]
    beta <- x[3:length(x)]
    
    tau <- diag((1-r)*(sigsq_e)/(sigsq_e+n*sigsq_v))
    D <- diag(sigsq_e/n)
    
    G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
    H <- G%*%X_bar-X_pop
    Hinv <- solve(t(H)%*%H)
    
    P <- diag(length(Y_bar))-H%*%Hinv%*%t(H)
    
    A<-(G-diag(length(Y_bar)))%*%Y_bar-H%*%beta
    
    
    return(sum(diag(D-2*tau%*%D))+t(A)%*%A)
  }
  
  
  result <- optim(par = c(sigsq_v,sigsq_e,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,n=n,r=r, lower = c(sigsq_v+0.00001,sigsq_e+0.00001,beta-0.15),upper = c(sigsq_v+0.00002,sigsq_e+0.00002,beta+0.15),method ="L-BFGS-B")
  
  sigsq_v <- result$par[1]
  sigsq_e <- result$par[2]
  beta <- result$par[3:length(result$par)]
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  
  theta <- X_pop%*%beta+G%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,beta=beta,theta=theta))
}

OBP_design <- function (X, X_pop, Y, r,sigsq_e, sigsq_v,beta){
  
  X_bar <- apply(X,3,rowMeans)
  #X_bar <- X_pop
  
  Y_bar <- rowMeans(Y)  
  n <- apply(X, 1, length)/dim(X_bar)[2]
  
  Qr <- function(x,X_bar,Y_bar,Y, X_pop,n,r){
    sigsq_e <- x[1]
    sigsq_v <- x[2]
    beta <- x[3:length(x)]
    N <-n/r
    
    ui.sq <- rowMeans(Y^2)-((apply(Y,1,sd))^2)*(N-1)/N
    G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
    H <- (X_pop-G%*%X_bar)
    L <- G^2%*%X_bar + (diag(length(Y_bar))-2*G)%*%X_pop
    Hinv <- solve(t(H)%*%H)
    
    P <- G^2-L%*%Hinv%*%t(L)
    #return(sum((diag(length(Y_bar))-2*G)%*%ui.sq)+(t(Y_bar)%*%P%*%Y_bar))
    return(t(beta)%*%t(X_pop-G%*%X_bar)%*%(X_pop-G%*%X_bar)%*%beta-2*t(Y_bar)%*%L%*%beta+t(Y_bar)%*%G^2%*%Y_bar+sum((diag(length(Y_bar))-2*G)%*%ui.sq))
  }
  
  
  result <- optim(par = c(sigsq_e,sigsq_v,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,Y=Y,n=n,r=r, lower = c(sigsq_e+0.00001,sigsq_v+0.00001,beta-0.15),upper = c(sigsq_e+0.00002,sigsq_v+0.00002,beta+0.15),method ="L-BFGS-B")
  
  sigsq_e <- result$par[1]
  sigsq_v <- result$par[2]
  beta <- result$par[3:length(result$par)]
  
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  H <- (X_pop-G%*%X_bar)
  L <- G^2%*%X_bar + (diag(length(Y_bar))-2*G)%*%X_pop
  Hinv <- solve(t(H)%*%H)
  
  
  #beta <- Hinv%*%t(L)%*%Y_bar
  theta <- X_pop%*%beta+diag(r+(1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v))%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
}


mle_est <- function (X, X_pop, Y, r){
  
  
  Y_bar <- rowMeans(Y)
  X_bar <- apply(X,3,rowMeans)
  n <- apply(X, 1, length)/dim(X_bar)[2]
  
  id <- seq(1:length(Y_bar))
  X_reshape<-c()
  for ( i in id){
    X_samp <- unlist(X[i,,])
    X_new <- cbind(rep(i,dim(X_samp)[1]),X_samp)
    X_reshape <- rbind(X_reshape,X_new)
  }
  Y_reshape <- c(t(Y))
  data_samp <- data.frame(Y=Y_reshape,id=X_reshape[,1],X_new1=X_reshape[,2],X_new2=X_reshape[,3])
  res <- lmer(Y~X_new1+X_new2+(1|id)-1,data=data_samp,REML = FALSE)
  v<-as.data.frame(VarCorr(res),comp="Variance")
  sigsq_v<-v$vcov[1]
  sigsq_e <- v$vcov[2]
  
  beta<- as.vector(summary(res)$coefficients[,1])
  

  theta <- X_pop%*%beta+diag(r+(1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v))%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
}

library(abind)
sim_data <- function(m,n,p,r,sigsq_e,sigsq_v){
  X <- array(rlnorm(m*n*p/r,1,0.5),c(m,n/r,p))
  beta = rep(0,p)
  v <- matrix(rep(rnorm(m,sd=sqrt(sigsq_v)),n/r),m,n/r)
  e <- matrix(rnorm(m*n/r,sd=sqrt(sigsq_e)),m,n/r)
  Y <- apply(X, 2, FUN = function(x){x%*%beta})+v+e+10
  #Y <- 10+v+e
  r <- rep(r,m)
  
  
  # X_samp <- array(NA,c(m,n,p))
  # Y_samp <- matrix(NA,m,n)
  # for (i in 1:m){
  #   index <- order(sqrt(X[i,,1]^2+X[i,,2]^2))[1:n]
  #   X_samp[i,,] <- X[i,index,]
  #   Y_samp[i,] <- Y[i,index]
  # }
  
  index <- sample(1:(n/r[1]),n)
  X_samp <- X[,index,]
  Y_samp <- Y[,index]
  return(list(X=X,v=v,e=e,Y=Y,X_samp=X_samp,Y_samp=Y_samp,r=r))
}

mse_fh_al  <- c()
mse_fh_uc  <- c()
mse_fh_model  <- c()
mse_ds <- c()
mse_obp <- c()
mse_dir <- c()
mse_mle <- c()
for (i in 1:100){
  data <-sim_data(m=40,n=4,p=2,r=1/250,sigsq_e = 6,sigsq_v = 1)
  
  sigsq_e_est_direct <- mean((apply(data$Y_samp,1,sd))^2)
  X_pop <- apply(data$X,3,rowMeans)
  
  res_mle <- mle_est(X=data$X_samp,X_pop=X_pop,Y=data$Y_samp,r=data$r)
  sigsq_e_est <-res_mle$sigsq_e
  sigsq_v_est <-res_mle$sigsq_v
  res_fh <- OBP_fh(X=data$X_samp,X_pop=X_pop,Y=data$Y_samp,r=data$r,sigsq_e = sigsq_e_est_direct,sigsq_v = res_mle$sigsq_v)
  res_fh_uc <- OBP_fh_uc(X=data$X_samp,X_pop=X_pop,Y=data$Y_samp,r=data$r,sigsq_e = sigsq_e_est_direct,sigsq_v = res_mle$sigsq_v)
  res_fh_model <- OBP_fh_model(X=data$X_samp,X_pop=X_pop,Y=data$Y_samp,r=data$r,sigsq_e = sigsq_e_est_direct,sigsq_v = res_mle$sigsq_v)
  
  res_new <- OBP_model(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta)
  res_ds <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta)

  
  
  theta_true <- apply(data$X, 2, FUN = function(x){x%*%rep(0,2)})+data$v[,1]+10
  #theta_true <- 10+data$v[,1]
  theta_true <- rowMeans(theta_true)
  Y_bar <- rowMeans(data$Y_samp)
  
  mse_fh_al[i] <- mean((theta_true-res_fh$theta)^2)
  mse_fh_uc[i] <- mean((theta_true-res_fh_uc$theta)^2)
  mse_fh_model[i] <- mean((theta_true-res_fh_model$theta)^2)
  mse_ds[i] <- mean((theta_true-res_ds$theta)^2)
  mse_obp[i] <- mean((theta_true-res_new$theta)^2)
  mse_dir[i] <- mean((theta_true-Y_bar)^2)
  mse_mle[i] <- mean((theta_true-res_mle$theta)^2)
}


mse <- data.frame(Type=c(rep('dir',100),rep('eblup',100),rep('fh_al',100),rep('fh_uc',100), rep('obp-ds',100),rep('obp',100),rep('fh-obp',100)),MSE=c(mse_dir,mse_mle,mse_fh_al,mse_fh_uc,mse_ds,mse_obp,mse_fh_model))
ggplot(mse, aes(x=Type, y=abs(MSE), fill=Type)) +geom_boxplot(outliers = FALSE)
saveRDS(mse,file ='C:/Users/jiang/Desktop/paper/sae/SAE_RES_new/mse_m100_tot10.rds' )
