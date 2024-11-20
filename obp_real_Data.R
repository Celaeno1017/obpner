BHF <- read.csv('D:/download/BHF_data.csv')
BHF[2,2] <-1
BHF[3,2] <-1
BHF$County[4:37] <-BHF$County[4:37]-2 
PopCorn_1 <- (BHF$PopCornPix[1]*BHF$Ncounty[1]+BHF$PopCornPix[2]*BHF$Ncounty[2]+BHF$PopCornPix[3]*BHF$Ncounty[3])/sum(BHF$Ncounty[1:3])
PopBean_1 <- (BHF$PopBeansPix[1]*BHF$Ncounty[1]+BHF$PopBeansPix[2]*BHF$Ncounty[2]+BHF$PopBeansPix[3]*BHF$Ncounty[3])/sum(BHF$Ncounty[1:3])
BHF[1,8] <- PopCorn_1
BHF[1,9] <- PopBean_1

BHF[1,3] <- BHF[1,3]+BHF[2,3]+BHF[3,3]

BHF[2,3] <- NA
BHF[3,3] <- NA
BHF[2,8] <- NA
BHF[3,8] <- NA
BHF[2,9] <- NA
BHF[3,9] <- NA

r <- table(BHF$County)/BHF$Ncounty[!is.na(BHF$Ncounty)]
Y_samp <- matrix(NA,10,6)
X_samp <- array(NA,c(10,6,3))
n<- rep(0,10)
for (i in 1:10){
  n[i] <- length(which(BHF$County==i))
  Y_samp[i,1:n[i]] <- BHF$CornHec[which(BHF$County==i)]
  
  X_samp[i,1:n[i],1] <- rep(1,length(BHF$SoyBeansPix[which(BHF$County==i)]))
  X_samp[i,1:n[i],2] <- BHF$CornPix[which(BHF$County==i)]
  X_samp[i,1:n[i],3] <- BHF$SoyBeansPix[which(BHF$County==i)]
 
  
}
X_pop <- BHF[!is.na(BHF$Ncounty),8:9]
X_pop <- cbind(rep(1,dim(X_pop)[1]),X_pop)

######Analysis


#####fh area-level model##############
OBP_fh <- function (X, X_pop,Y, r,n,sigsq_e, sigsq_v){
  
  X_bar <- X_pop
  #X_bar <- apply(X,3,rowMeans)
  Y_bar <- rowMeans(Y,na.rm=TRUE)  
  
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
OBP_fh_uc <- function (X, X_pop,Y, r,n, sigsq_e, sigsq_v){
  
  X_bar <- X_pop
  X_unit <- X_bar[rep(row.names(X_bar), times = n),]
  n_ori <- n
  n <- rep(1,sum(n))
  
  Y_unit <- as.vector(t(Y))
  Y_unit <- Y_unit[!is.na(Y_unit)]
  
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
  theta <- X_pop%*%beta 
  add<- diag((sigsq_v)/(sigsq_e/n+sigsq_v))%*%(Y_unit-(X_unit)%*%beta)
  
  splitAt <- function(x, pos) {
    unname(split(x, cumsum(seq_along(x) %in% pos)))
  }
  
  theta<- theta+sapply(splitAt(add, c(1, cumsum(n_ori)+1)), mean)
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
  
  return(result)
}


OBP_fh_model <- function (X, X_pop,Y, r,n,sigsq_e, sigsq_v){
  
  X_bar <- X_pop
  #X_bar <- apply(X,3,rowMeans)
  Y_bar <- rowMeans(Y,na.rm=TRUE)  
  
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
  
  X_bar <- apply(X,3,rowMeans,na.rm=TRUE)
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  H <- G%*%X_bar-X_pop
  Hinv <- solve(t(H)%*%H)
  
  
  beta <- Hinv%*%t(H)%*%(G-diag(length(Y_bar)))%*%Y_bar
  
  theta <- X_pop%*%beta+G%*%(Y_bar-X_bar%*%beta)
  
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
  
  return(result)
}



OBP_est_new <- function (X, X_pop, Y, r,n,sigsq_e, sigsq_v){
  X_bar <- apply(X,3,rowMeans,na.rm=TRUE)
  Y_bar <- rowMeans(Y,na.rm=TRUE)  
  
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

OBP_model <- function (X, X_pop, Y, r,n,sigsq_e, sigsq_v,beta){
  X_bar <- apply(X,3,rowMeans,na.rm=TRUE)
  Y_bar <- rowMeans(Y,na.rm=TRUE)  
  
  Qr <- function(x,X_bar,Y_bar,X_pop,sigsq_e,n,r){
    sigsq_v <- x[1]
    #sigsq_e <- x[2]
    beta <- x[2:length(x)]
    
    tau <- diag((1-r)*(sigsq_e)/(sigsq_e+n*sigsq_v))
    D <- diag(sigsq_e/n)
    
    G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
    H <- G%*%X_bar-X_pop
    Hinv <- solve(t(H)%*%H)
    
    P <- diag(length(Y_bar))-H%*%Hinv%*%t(H)
    
    A<-(G-diag(length(Y_bar)))%*%Y_bar-H%*%beta
    
    
    return(sum(diag(D-2*tau%*%D))+t(A)%*%A)
  }
  
  
  result <- optim(par = c(sigsq_v,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,sigsq_e=sigsq_e,n=n,r=r, lower = c(sigsq_v*0.95,beta-abs(beta)*0.5),upper = c(sigsq_v*1.05,beta+abs(beta)*0.5),method ="L-BFGS-B")
  
  sigsq_v <- result$par[1]
  #sigsq_e <- result$par[2]
  beta <- result$par[2:length(result$par)]
  G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
  
  theta <- X_pop%*%beta+G%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,beta=beta,theta=theta))
}

OBP_design <- function (X, X_pop, Y, r,n,sigsq_e, sigsq_v,beta){
  
  X_bar <- apply(X,3,rowMeans,na.rm=TRUE)
  #X_bar <- X_pop
  
  Y_bar <- rowMeans(Y,na.rm=TRUE)  
  
  Qr <- function(x,X_bar,Y_bar,Y, X_pop,n,r){
    sigsq_e <- x[1]
    sigsq_v <- x[2]
    beta <- x[3:length(x)]
    N <-n/r
    
    ui.sq <- rowMeans(Y^2,na.rm=TRUE)-((apply(Y,1,sd,na.rm=TRUE))^2)*(N-1)/N
    G <- diag((1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v)+r)
    H <- (X_pop-G%*%X_bar)
    L <- G^2%*%X_bar + (diag(length(Y_bar))-2*G)%*%X_pop
    Hinv <- solve(t(H)%*%H)
    
    P <- G^2-L%*%Hinv%*%t(L)
    #return(sum((diag(length(Y_bar))-2*G)%*%ui.sq)+(t(Y_bar)%*%P%*%Y_bar))
    return(t(beta)%*%t(X_pop-G%*%X_bar)%*%(X_pop-G%*%X_bar)%*%beta-2*t(Y_bar)%*%L%*%beta+t(Y_bar)%*%G^2%*%Y_bar+sum((diag(length(Y_bar))-2*G)%*%ui.sq))
  }
  
  
  result <- optim(par = c(sigsq_e,sigsq_v,beta), fn = Qr, X_bar = X_bar,Y_bar=Y_bar,X_pop=X_pop,Y=Y,n=n,r=r, lower = c(0.95*sigsq_e,0.95*sigsq_v,beta-abs(beta)*0.1),upper = c(1.01*sigsq_e,1.01*sigsq_v,beta+abs(beta)*0.1),method ="L-BFGS-B")
  
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


mle_est <- function (X, X_pop, Y, r,n){
  
  
  Y_bar <- rowMeans(Y,na.rm=TRUE)
  X_bar <- apply(X,3,rowMeans,na.rm=TRUE)
  n <- n
  
  id <- seq(1:length(Y_bar))
  X_reshape<-c()
  for ( i in id){
    X_samp <- unlist(X[i,,])
    X_new <- cbind(rep(i,dim(X_samp)[1]),X_samp)
    X_reshape <- rbind(X_reshape,X_new)
  }
  Y_reshape <- c(t(Y))
  data_samp <- data.frame(Y=Y_reshape,id=X_reshape[,1],X_new1=X_reshape[,2],X_new2=X_reshape[,3],X_new3=X_reshape[,4])
  res <- lmer(Y~X_new1+X_new2+X_new3+(1|id)-1,data=data_samp,REML = FALSE)
  v<-as.data.frame(VarCorr(res),comp="Variance")
  sigsq_v<-v$vcov[1]
  sigsq_e <- v$vcov[2]
  
  beta<- as.vector(summary(res)$coefficients[,1])
  
  
  theta <- X_pop%*%beta+diag(r+(1-r)*n*sigsq_v/(sigsq_e+n*sigsq_v))%*%(Y_bar-X_bar%*%beta)
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
}


mle_est_uc <- function (X, X_pop, Y, r,n){
  X_bar <- X_pop
  
  Y_bar <- rowMeans(Y,na.rm=TRUE)
  id <- seq(1:length(Y_bar))
  
  X_reshape <- X_bar[rep(row.names(X_bar), times = n),]
  X_reshape <- cbind(rep(id,times=n),X_reshape)
  
  
  
  
 
  
  Y_reshape <- c(t(Y))
  Y_reshape <- Y_reshape[!is.na(Y_reshape)]
  data_samp <- data.frame(Y=Y_reshape,id=X_reshape[,1],X_new1=X_reshape[,2],X_new2=X_reshape[,3],X_new3=X_reshape[,4])
  res <- lmer(Y~X_new1+X_new2+X_new3+(1|id)-1,data=data_samp)
  v<-as.data.frame(VarCorr(res),comp="Variance")
  sigsq_v<-v$vcov[1]
  sigsq_e <- v$vcov[2]
  
  beta<- as.vector(summary(res)$coefficients[,1])
  
  Y_unit <- Y_reshape
  X_unit <- X_bar[rep(row.names(X_bar), times = n),]
  
  
  n <- rep(1,sum(n))
  
  theta <- X_pop%*%beta + .colMeans(diag((sigsq_v)/(sigsq_e/n+sigsq_v))%*%(Y_unit-(X_unit)%*%beta),length(n)/length(apply(X, 1, length)), length(apply(X, 1, length)))
  
  
  return(list(sigsq_v=sigsq_v,sigsq_e=sigsq_e,beta=beta,theta=theta))
}




library(lme4)



sigsq_e_est_direct <- mean((apply(Y_samp,1,sd,na.rm=TRUE))^2)
X_pop <- as.matrix(X_pop)
res_mle <- mle_est(X=X_samp,X_pop=X_pop,Y=Y_samp,r=r,n=n)
res_mle_uc <- mle_est_uc(X=X_samp,X_pop=X_pop,Y=Y_samp,r=r,n=n)

sigsq_e_est <-res_mle$sigsq_e
sigsq_v_est <-res_mle$sigsq_v
res_fh <- OBP_fh(X=X_samp,X_pop=X_pop,Y=Y_samp,r=r,n=n,sigsq_e = sigsq_e_est_direct,sigsq_v = res_mle$sigsq_v)
res_fh_uc <- OBP_fh_uc(X=X_samp,X_pop=X_pop,Y=Y_samp,r=r,n=n,sigsq_e = sigsq_e_est_direct,sigsq_v = res_mle$sigsq_v)
res_fh_model <- OBP_fh_model(X=X_samp,X_pop=X_pop,Y=Y_samp,r=r,n=n,sigsq_e = sigsq_e_est_direct,sigsq_v = res_mle$sigsq_v)

res_new <- OBP_model(X=X_samp,Y=Y_samp,X_pop=X_pop,r=r,n=n,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta)
res_ds <- OBP_design(X=X_samp,Y=Y_samp,X_pop=X_pop,r=r,n=n,sigsq_e = sigsq_e_est,sigsq_v = res_mle$sigsq_v,beta=res_mle$beta)


