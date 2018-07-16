# https://docs.google.com/document/d/1YHD1ZskYSWQ60HCj_OChLxomAuNcERHPhEOeLs0qSMc/edit
# Mandar os Dados  rafaerbisti@gmail.com
# E-mail Prof Helio Migon
# migon@im.ufrj.br

#gustavolibel@gmail.com

# Códigos da Semana 2 -------
# Autora: Raíra Marotta
# e-mail: raira@dme.ufrj.br
# Estes códigos estão disponibilizados em https://github.com/rairamarotta


# Modelo Normal com variância desconhecida ----------


normal_var_unkw <- function(y, m0, C0, n0, d0, W, F, G, D){
  
  #Auxilio Desconto
  if(is.null(dim(G)) == FALSE){
    matrixaux <- G == 0 
    if(G[1,1] == 1 & G[1,2] == 1 & G[2,2] == 1 ) {matrixaux[2] <- 0}
    tira <- which(matrixaux == 1)
    mantem <- which(matrixaux == 0)}
  
  #Definindo os objetos
  n <-nrow(F) ; r <-1
  T <- length(y)
  S0 <- d0/n0
  mt <- matrix(0,nrow=n,ncol=T+1)
  Ct <- array(rep(diag(n),T+1),dim=c(n,n,T+1))
  Rt <- array(rep(0,T+1),dim=c(n,n,T+1))
  Pt <- array(rep(0,T+1),dim=c(n,n,T+1))
  Wt <- array(rep(0,T+1),dim=c(n,n,T+1))
  ft <- matrix(0,nrow=T+1,ncol=r)
  at <- matrix(0,nrow=n,ncol=T+1)
  Qt <- matrix(0,nrow=T+1,ncol=r)
  et <- matrix(0,nrow=T+1,ncol=r)
  At <- matrix(0,nrow=n,ncol=T+1)
  dt <- matrix(0,nrow=T+1,ncol=r)
  nt <- matrix(0,nrow=T+1,ncol=r)
  St <- matrix(0,nrow=T+1,ncol=r)
  LSt <- matrix(0,nrow=T+1,ncol=r)
  LIt <- matrix(0,nrow=T+1,ncol=r)
  ## Passo t=1
  
  if(is.null(D) == TRUE){
    at[,1] <- G%*%m0
    Rt[,,1] <- G%*%C0%*%(t(G)) + W} else{
      if(is.null(dim(D)) == TRUE){
        at[,1] <- G%*%m0
        Pt[,,1] <- G%*%C0%*%(t(G))
        Wt[,,1] <- (D^2 -1)*Pt[,,1]
        Rt[,,1] <- D%*%Pt[,,1]%*%t(D)
      } else{
        at[,1] <- G%*%m0
        Pt[,,1] <- G%*%C0%*%(t(G))
        aux1 <- (D^2)*Pt[,,1]
        aux2 <- Pt[,,1]
        aux.W <- D^2- matrix(1,ncol=n, nrow = n)
        
        Wt[,,1] <-  Pt[,,1]*aux.W
        for(i in 1: length(mantem)){
          aux2[mantem[i]] <- 0}
        Rt[,,1] <- aux1 + aux2}}
  
  # Previsăo 1 passo-a-frente
  
  ft[1,] <- t(F[,1])%*%at[,1]
  Qt[1,] <- t(F[,1])%*%Rt[,,1]%*%F[,1]+S0
  
  LSt[1] <- ft[1,] + qt(0.975, n0)*sqrt(Qt[1,]) 
  LIt[1] <- ft[1,] + qt(0.025, n0)*sqrt(Qt[1,]) 
  
  # Posteriori em t=1
  
  At[,1] <- Rt[,,1]%*%F[,1]*(1/Qt[1,])
  et[1,] <- y[1]-ft[1,]
  
  nt[1] <- n0+1
  dt[1] <- d0+S0*(et[1,]^2)/Qt[1,]
  St[1] <- dt[1]/nt[1]
  
  mt[,1] <- at[,1]+At[,1]*et[1,]
  Ct[,,1] <- (St[1]/S0)*(Rt[,,1]-At[,1]%*%t(At[,1])*Qt[1,])
  
  
  for(t in 2:(T+1)){			
    
    # Priori em t
    
    if(is.null(D) == TRUE){
      at[,t] <- G%*%mt[,t-1]
      Rt[,,t] <- G%*%Ct[,,t-1]%*%(t(G)) + W} else{
        if(is.null(dim(D)) == TRUE){
          at[,t] <- G%*%mt[,t-1]
          Pt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))
          Wt[,,t] <- (D^2 -1)*Pt[,,t]
          Rt[,,t] <- D%*%Pt[,,t]%*%t(D)
        } else{
          at[,t] <- G%*%mt[,t-1]
          Pt[,,t] <- G%*%Ct[,,t-1]%*%(t(G))
          aux1 <- (D^2)*Pt[,,t]
          aux2 <- Pt[,,t]
          Wt[,,t] <-  Pt[,,t]*aux.W
          for(i in 1: length(mantem)){
            aux2[mantem[i]] <- 0}
          Rt[,,t] <- aux1 + aux2}}
    
    
    # Previsăo 1 passo-a-frente
    ft[t,] <- t(F[,t])%*%at[,t]
    Qt[t,] <- t(F[,t])%*%Rt[,,t]%*%F[,t]+St[t-1]
    
    LSt[t] <- ft[t,] + qt(0.975, nt[t-1])*sqrt(Qt[t,]) 
    LIt[t] <- ft[t,] + qt(0.025, nt[t-1])*sqrt(Qt[t,]) 
    
    # Posteriori em t
    At[,t] <- Rt[,,t]%*%F[,t]*(1/Qt[t,])
    et[t,] <- y[t]-ft[t,]
    
    nt[t] <- nt[t-1]+1
    dt[t] <- dt[t-1]+St[t-1]*(et[t,]^2)/Qt[t,]
    St[t] <- dt[t]/nt[t]
    
    mt[,t] <- at[,t]+At[,t]*et[t,]
    Ct[,,t] <- (St[t]/St[t-1])*(Rt[,,t]-At[,t]%*%t(At[,t])*Qt[t,])
  }
  
  if(is.null(W) == TRUE){
    result <- list(mt,Ct,ft,Qt, et, Pt, Wt, LSt, LIt, nt, Rt, at, St)
    names(result) <- c("mt", "Ct", "ft", "Qt", "et", "Pt", "Wt", "LSt", "LIt", "nt", "Rt", "at","St")} else{
      result <- list(mt,Ct,ft,Qt, et, LSt, LIt, nt, Rt, at,St)  
      names(result) <- c("mt", "Ct", "ft", "Qt", "et", "LSt", "LIt", "nt", "Rt", "at","St")
    }
  return(result)
}

library(dlm)

#library(readxl)
#setwd("D://3aEscMat")
df_usp <- read.csv(file="USP.csv", sep=",")
y <- df_usp$freq

m0 <- 0; C0<- 1000; W<- 40; n0 <- 1; d0 <- 1

# Definindo a matriz F
F1 <- matrix(1,nrow=1,ncol=(length(y)+1))
F <- F1
F

resultados <-normal_var_unkw(y, m0, C0, n0, d0, W, F, G = 1, D = NULL)

# Previsão a 1 passo
plot(y,pch = 20, ylab = " ", xlab = " ")
lines(resultados$ft, col = 2)
lines(resultados$LSt, col = 4, lty = 2)
lines(resultados$LIt, col = 4, lty = 2)
lines(c(0,rep(c(rep(c(2000),12),rep(c(0),12)),1)), col = 3, lty = 2)
legend("topright", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

#y <- as.numeric(AirPassengers)

m0 <- 0; C0<- 1000; D <- 1/sqrt(0.6) ; n0 <- 1; d0 <- 1

# Definindo a matriz F
F1 <- matrix(1,nrow=1,ncol=(length(y)+1))
F <- F1
F

resultados <-normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G = 1, D = D)

# Previsão a 1 passo
plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
lines(resultados$LSt, col = 4, lty = 2)
lines(resultados$LIt, col = 4, lty = 2)
grid.points(x=12)
legend("topright", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Modelo de Segunda  ordem ---------

# Definindo a matriz F
n <- 2
F1 <- matrix(c(1,0),nrow=n,ncol=(length(y)+1))
F <- F1
F

# Definindo a matriz G
G <- matrix(c(1,0,1,1),2,2)
G

# Definindo as quantidades iniciais
m0 <- c(0,0)
C0 <- diag(100,n,n) 
W = diag(c(0.02, 0.01))

resultados <- normal_var_unkw(y, m0, C0, n0, d0, W, F, G, D = NULL)

# Estimação do nível

LS <- resultados$mt[1,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[1,1,])
LI <- resultados$mt[1,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[1,1,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(min(LI) - 1, max(LS)+1))
lines(resultados$mt[1,], col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Previsão a 1 passo
LS <- dropFirst(qnorm(0.975, resultados$ft, sqrt(resultados$Qt)))
LI <- dropFirst(qnorm(0.025, resultados$ft, sqrt(resultados$Qt)))

plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Fator de crescimento
LS <- resultados$mt[2,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[2,2,])
LI <- resultados$mt[2,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[2,2,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[2,], type = "l", col =2, ylim= c(-2,2))
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)


# Descontos -------------

# Desconto 1
m0 <- 0; C0<- 100; W = NULL
F <- matrix(1,nrow=1,ncol=(length(y)+1))
G <- 1
D <- 1
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsão a 1 passo
plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Desconto 0.9
m0 <- 0; C0<- 1000; W = NULL
F <- matrix(1,nrow=1,ncol=(length(y)+1))
G <- 1
D <- 1/sqrt(0.9)
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsão a 1 passo

plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Desconto 0.8
m0 <- 0; C0<- 1000; W = NULL
F <- matrix(1,nrow=1,ncol=(length(y)+1))
G <- 1
D <- 1/sqrt(0.8)
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsão a 1 passo

plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Regressão múltipla ---------


# Modelo
n <- 2
m0  <- c(0,0); C0<- diag(100,2); W = NULL
T <- length(y)

#Definindo F
F <- matrix(c(rep(1, length(y)),df_usp$freq),ncol=T, nrow=n, byrow=T)
F <- cbind(F, c(1,1)) #gambiarra no código (pq antes estava indo até t+1)

#Definindo G
G <- diag(1,2)

#Definindo D
D1 <- 1/sqrt(0.9)
D2 <- 1/sqrt(0.99)
D<- bdiag(D1,D2)

# Resultados
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsao
plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(3,15))
lines(dropFirst(resultados$ft), col = 2)
lines(resultados$LSt, col = 4, lty = 2)
lines(resultados$LIt, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Efeito do Preço

LS <- resultados$mt[2,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[2,2,])
LI <- resultados$mt[2,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[2,2,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

par(mfrow=c(2,1), mar = c(4,2,2,1))
plot(resultados$mt[2,], type = "l", col =2, ylim=c(-5,3),
     xlab = expression(paste("E[",beta[t],"|D"[t], "]")))
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Efeito do preço x Preço
plot(na.omit(resultados$mt[2,])*df_usp$freq, type = "l", col =2, ylim=c(-5,3),
     xlab = expression(paste("E[",beta[t],"|D"[t], "] . x"[t])))
lines((df_usp$freq)*LS, col = 4, lty = 2)
lines((df_usp$freq)*LI, col = 4, lty = 2)



# Modelo com 1 harmonico
# Definindo as quantidades
n <- 4
T <- length(y)
w <- (2*pi)/(12)
G1 <- diag(1,2)
G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
G <- dlm::bdiag(G1, G2)

F <- matrix(c(rep(1, length(y)),df_usp$freq, 
              rep(1, length(y)),  rep(0, length(y))),ncol=T, nrow=n, byrow=T)
F <- cbind(F, c(1,1,1,1)) #gambiarra no código (pq antes estava indo até t+1)

#Desconto da tendência
D1 <- 1/sqrt(0.97)

#Desconto da covariável
D2 <- 1/sqrt(0.99)

#Desconto da sazonalidade
D3 <- matrix(1/sqrt(0.99),2,2)
D <- dlm::bdiag(D1,D2, D3)

F.prev <- F
F.prev[2,] <- 0 #Não temos mais dados da covariável, por isso, zero
F.prev <- F.prev[,1:12]
G.prev <- G

m0 <- rep(0,n)
C0 <- diag(1000,n)

#Resultados
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsao
plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
lines(resultados$LSt, col = 4, lty = 2)
lines(resultados$LIt, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Efeito do harmonico

LS <- resultados$mt[3,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[3,3,])
LI <- resultados$mt[3,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[3,3,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[3,], type = "l", col =2, xlab = " ", ylab = " ")
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)


# Modelo com 1 harmonico NOVO

# Definindo as quantidades
n <- 4
w <- (2*pi)/(12)
G1 <- matrix(c(1,1,0,1),ncol=2,nrow=2,byrow=TRUE)
G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
G <- bdiag(G1, G2)
T <- length(y)

F <- matrix(c(rep(1,length(y)),rep(0,length(y)),rep(1,length(y)),rep(0,length(y))),ncol=T,nrow=n,byrow=TRUE)
F <- cbind(F,c(1,0,1,0)) #gambiarra no código (pq antes estava indo até t+1)

#Desconto da tendência
D1 <- matrix(1/sqrt(0.8),2,2)

#Desconto da sazonalidade
D3 <- matrix(1/sqrt(0.95),2,2)
D <- bdiag(D1,D3)

F.prev <- F
G.prev <- G

m0 <- rep(0,n)
C0 <- diag(100,n)

n0 <- 10
d0 <- 1

#Resultados
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsao
plot(y,pch = 20, ylab = "Porcentagem de acesso", xlab = "Meses", ylim= c(0,110),main="Modelo com 1 Harmônico")
lines(dropFirst(resultados$ft), col = 2)
lines(resultados$LSt, col = 4, lty = 2)
lines(resultados$LIt, col = 4, lty = 2)
legend("topright", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)


# Efeito do harmonico

LS <- resultados$mt[3,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[3,3,])
LI <- resultados$mt[3,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[3,3,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[3,], type = "l", col =2, xlab = " ", ylab = " ",main="Efeito do Harmônico")
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Modelo com 2 harmonicos ------


# Definindo as quantidades
n <- 6
w <- (2*pi)/(12)
G2 <- matrix(c(cos(w),-sin(w),sin(w),cos(w)),2,2)
G3 <- matrix(c(cos(2*w),-sin(2*w),sin(2*w),cos(2*w)),2,2)
G <- dlm::bdiag(G1, G2, G3)


F <- matrix(c(rep(1, length(y)),df_usp$freq, 
              rep(1, length(y)),  rep(0, length(y)),
              rep(1, length(y)),  rep(0, length(y))),ncol=T, nrow=n, byrow=T)
F <- cbind(F, c(1,1,1,1,1,1)) #gambiarra no código (pq antes estava indo até t+1)

#Desconto da tendência
D1 <- 1/sqrt(0.97)

#Desconto da covariável
D2 <- 1/sqrt(0.99)

#Desconto da sazonalidade
D3 <- matrix(1/sqrt(0.99),2,2)
D4 <- matrix(1/sqrt(0.99),2,2)

#Desconto final
D <- dlm::bdiag(D1,D2,D3,D4)

#Previsão
F.prev <- F
#F.prev[2,] <- 0 #Não temos mais dados da covariável, por isso, zero
#F.prev <- F.prev[,1:12]
G.prev <- G

m0 <- rep(0,n)
C0 <- diag(100,n)

#Resultados
resultados <- normal_var_unkw(y, m0, C0, n0, d0, W = NULL, F, G, D)

# Previsao
plot(y,pch = 20, ylab = " ", xlab = " ")
lines(dropFirst(resultados$ft), col = 2)
lines(resultados$LSt, col = 4, lty = 2)
lines(resultados$LIt, col = 4, lty = 2)
legend("topleft", legend = c(expression(y[t]),expression(f[t])),
       bty = "n", pch = c(20, NA), lty = c(NA,1), col = c(1,2),
       cex = 1.2)

# Efeito do harmonico

LS <- resultados$mt[3,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[3,3,])
LI <- resultados$mt[3,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[3,3,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

plot(resultados$mt[3,], type = "l", col =2, xlab = " ", ylab = " ")
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)

# Análise Retrospectiva ---------
analise.restrospectiva <- function(mt, Ct, Rt, at){
  
  T <- ncol(mt) -1
  ms1 <- matrix(0,nrow=n,ncol=T)
  Rs1 <- array(rep(diag(n),T),dim=c(n,n,T))
  
  ms1[,T] <- mt[,T]
  Rs1[,,T] <- Ct[,,T]
  
  for(t in (T-1):1){
    ms1[,t] <- mt[,t] + Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(ms1[,t+1] - at[,t+1])
    Rs1[,,t] <- Ct[,,t] + Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1])%*%(Rt[,,t+1] - Rs1[,,t+1])%*%t(Ct[,,t]%*%t(G)%*%solve(Rt[,,t+1]))
  }
  result <- list(ms1, Rs1)
  names(result) <- c("ms1", "Rs1")
  return(result)
}

ms1 <- matrix(0,nrow=n,ncol=T)
Rs1 <- array(rep(diag(n),T),dim=c(n,n,T))

ms1[,T] <- resultados$mt[,T]
Rs1[,,T] <- resultados$Ct[,,T]

library(MASS)
for(t in (T-1):2){
  ms1[,t] <- resultados$mt[,t] + resultados$Ct[,,t]%*%G%*%ginv(resultados$Rt[,,t+1])%*%(ms1[,t+1] - resultados$at[,t+1])
  Rs1[,,t] <- resultados$Ct[,,t] + resultados$Ct[,,t]%*%G%*%ginv(resultados$Rt[,,t+1])%*%(Rs1[,,t+1]- resultados$Rt[,,t+1])%*%t(resultados$Ct[,,t]%*%G%*%(1/resultados$Rt[,,t+1]))
}


LS <- resultados$mt[1,] + qt(0.975, resultados$nt)*sqrt(resultados$Ct[1,1,])
LI <- resultados$mt[1,] + qt(0.025, resultados$nt)*sqrt(resultados$Ct[1,1,])
LS <- LS[1:length(LS)-1]
LI <- LI[1:length(LI)-1]

LSs <- ms1[1,] + qt(0.975, resultados$nt[-1])*sqrt(Rs1[1,1,])
LIs <- ms1[1,] + qt(0.025, resultados$nt[-1])*sqrt(Rs1[1,1,])
LSs <- LS[1:length(LS)-1]
LIs <- LI[1:length(LI)-1]



plot(y,pch = 20, ylab = " ", xlab = " ", ylim= c(min(LI) - 1, max(LS)+1))
lines(resultados$mt[1,], col = 2)
lines(ms1[1,], col = 4)
lines(LS, col = 2, lty = 2)
lines(LI, col = 2, lty = 2)

lines(LSs, col = 4, lty = 2)
lines(LIs, col = 4, lty = 2)