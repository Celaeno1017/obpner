###estimate OBP under NER model############
library(lme4)
####input X: 3-D array with i: area i; j: sample j in n_i; k: covariate k in 1-p.
####input Y: 2-D array with i: arear i; j; sample j in n_i.
####input r: sampling fraction for each area i


## constrained design based obp
OBP_design <- function (X, X_pop, Y, r,sigsq_e, sigsq_v,beta,delta){
  
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
  
  
  result <- optim(par = c(sigsq_e,sigsq_v,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,Y=Y,n=n,r=r, lower = c(sigsq_e+0.00001,sigsq_v+0.00001,beta*(1-delta)),upper = c(sigsq_e+0.00002,sigsq_v+0.00002,beta*(1+delta)),method ="L-BFGS-B")
  
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


## unconstrained design based obp
OBP_design_uncons <- function (X, X_pop, Y, r,sigsq_e, sigsq_v,beta){
  
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
  
  
  result <- optim(par = c(sigsq_e,sigsq_v,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,Y=Y,n=n,r=r,method ="L-BFGS-B")
  
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

## constrained model based obp
OBP_model <- function (X, X_pop, Y, r,sigsq_e, sigsq_v,beta,delta){
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
  
  
  result <- optim(par = c(sigsq_v,sigsq_e,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,n=n,r=r, lower = c(sigsq_v+0.00001,sigsq_e+0.00001,beta*(1-delta)),upper = c(sigsq_v+0.00002,sigsq_e+0.00002,beta*(1+delta)),method ="L-BFGS-B")
  
  sigsq_v <- result$par[1]
  sigsq_e <- result$par[2]
  beta <- result$par[3:length(result$par)]
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  
  theta <- X_pop%*%beta+G%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,beta=beta,theta=theta))
}

##eblup
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
  #e <- matrix(rnorm(m*n/r,sd=sqrt(sigsq_e)),m,n/r)
  sigsq_e <- rgamma(m,3,0.5)
  e <- t(sapply(sigsq_e,function(x){rnorm(n/r,sd=sqrt(x))}))
  Y <- apply(X, 2, FUN = function(x){x%*%beta})+v+e+5
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

mse_mle <- c()
mse_ds_un <- c()
mse_ds_0.1 <- c()
mse_ds_0.2 <- c()
mse_ds_0.3 <- c()
mse_ds_0.5 <- c()
mse_ds_0.7 <- c()
mse_ds_0.9 <- c()
mse_ds_1 <- c()
mse_ds_2 <- c()
mse_ds_5 <- c()
mse_ds_10 <- c()
mse_ds_0.05 <- c()
mse_ds_0.01 <- c()
mse_ds_0.005 <- c()

for (i in 1:100){
  data <-sim_data(m=400,n=4,p=2,r=1/250,sigsq_e = 6,sigsq_v = 1)
  
  sigsq_e_est_direct <- mean((apply(data$Y_samp,1,sd))^2)
  X_pop <- apply(data$X,3,rowMeans)
  
  res_mle <- mle_est(X=data$X_samp,X_pop=X_pop,Y=data$Y_samp,r=data$r)
  
  sigsq_e_est <-res_mle$sigsq_e
  sigsq_v_est <-res_mle$sigsq_v
  
  res_ds_0.1 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.1)
  res_ds_0.2 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.2)
  res_ds_0.3 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.3)
  res_ds_0.5 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.5)
  res_ds_0.7 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.7)
  res_ds_0.9 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.9)
  res_ds_1 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=1)
  res_ds_2 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=2)
  res_ds_5 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=5)
  res_ds_10 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=10)
  res_ds_0.05 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.05)
  res_ds_0.01 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.01)
  res_ds_0.005 <- OBP_design(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta,delta=0.005)
  res_ds_un <- OBP_design_uncons(X=data$X_samp,Y=data$Y_samp,X_pop=X_pop,r=data$r,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta)
  
  
  
  
  
  theta_true <- apply(data$X, 2, FUN = function(x){x%*%rep(0,2)})+data$v[,1]+5
  #theta_true <- 10+data$v[,1]
  theta_true <- rowMeans(theta_true)
  Y_bar <- rowMeans(data$Y_samp)
  
 
  mse_ds_0.1[i] <- mean((theta_true-res_ds_0.1$theta)^2)
  mse_ds_0.2[i] <- mean((theta_true-res_ds_0.2$theta)^2)
  mse_ds_0.3[i] <- mean((theta_true-res_ds_0.3$theta)^2)
  mse_ds_0.5[i] <- mean((theta_true-res_ds_0.5$theta)^2)
  mse_ds_0.7[i] <- mean((theta_true-res_ds_0.7$theta)^2)
  mse_ds_0.9[i] <- mean((theta_true-res_ds_0.9$theta)^2)
  mse_ds_1[i] <- mean((theta_true-res_ds_1$theta)^2)
  mse_ds_2[i] <- mean((theta_true-res_ds_2$theta)^2)
  mse_ds_5[i] <- mean((theta_true-res_ds_5$theta)^2)
  mse_ds_10[i] <- mean((theta_true-res_ds_10$theta)^2)
  mse_ds_0.05[i] <- mean((theta_true-res_ds_0.05$theta)^2)
  mse_ds_0.01[i] <- mean((theta_true-res_ds_0.01$theta)^2)
  mse_ds_0.005[i] <- mean((theta_true-res_ds_0.005$theta)^2)
  mse_ds_un[i] <- mean((theta_true-res_ds_un$theta)^2)
  mse_mle[i] <- mean((theta_true-res_mle$theta)^2)
  print(i)
}

mse <- data.frame(Type=c(rep('un',100),rep('10',100),rep('5',100),rep('2',100),rep('1',100),rep('0.9',100),rep('0.7',100),rep('0.5',100),rep('0.3',100),rep('0.2',100),rep('0.1',100),rep('0.05',100), rep('0.01',100),rep('0.005',100),rep('mle',100)),MSE=c(mse_ds_un,mse_ds_10,mse_ds_5,mse_ds_2,mse_ds_1,mse_ds_0.9,mse_ds_0.7,mse_ds_0.5,mse_ds_0.3,mse_ds_0.2,mse_ds_0.1,mse_ds_0.05,mse_ds_0.01,mse_ds_0.005,mse_mle))
mse$Type <- factor(mse$Type,levels = c('mle','0.005','0.01','0.05','0.1','0.2',
                                       '0.3','0.5','0.7','0.9',
                                       '1','2','5','10','un'))
levels(mse$Type)[levels(mse$Type)=="mle"] <- "EBLUP"
levels(mse$Type)[levels(mse$Type)=="un"] <- "unconstrained"
ggplot(mse, aes(x=Type, y=abs(MSE), fill=Type)) +geom_boxplot(outliers = TRUE)+
  xlab("Delta") + theme_bw()+theme(legend.position="none",axis.title.y=element_blank())
saveRDS(mse,file ='C:/Users/jiang/Desktop/paper/sae/SAE_RES_new/mse_m100_full5_comp_delta.rds' )

mse<-readRDS(file ='C:/Users/jiang/Desktop/paper/sae/SAE_RES_new/mse_m400_partial5_comp_delta.rds')


