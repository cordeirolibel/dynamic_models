library(dlm)
#library(readxl)

#Lendo CSV
df_usp <- read.csv(file="USP.csv", sep=",")

#Grafico de todo o dado
plot(df_usp$freq,type="l",ylab="freq(%)",xlab="time(mes)",main="USP (2004-)")

#Grafico ultimos 5 anos
df_usp5 <- tail(df_usp,n=12*5)
plot(df_usp5$freq,type="l",ylab="freq(%)",xlab="time(mes)",main="USP (Ultimos 5 anos)")


#===========================================
#=============== Modelos 
#===========================================

yt <- df_usp$freq

#---------------------------------
# ====> Modelo Dinamico de 1 Ordem

# Passo 1: Hyperparameters
m0 <- 50
C0 <- 10
W <- 0.01
V <- 0.1

# Passo 2: Primeira Interacao
# Evolucao
at <- m0
Rt <- C0+W
# Predicao
ft <- m0 #at ??
Qt <- C0+W+V #Rt+V ??
# Atualizacao
At <- Rt/Qt
mt <- m0+Rt*(yt-m0)/Qt
Ct <- Rt-At^2*Qt

# Passo 3: Demais Iteracoes
for(t in 2:length(yt)){
  # Evolucao
  at[t] <- mt[t-1]
  Rt[t] <- Ct[t-1]+W
  # Predicao
  ft[t] <- mt[t-1]
  Qt[t] <- Ct[t-1]+W+V
  # Atualizacao
  At[t] <- Rt[t]/Qt[t]
  mt[t] <- mt[t-1]+Rt[t]*(yt[t]-mt[t-1])/Qt[t]
  Ct[t] <- Rt[t]-At[t]^2*Qt[t] 
}

# Passo 4: Grafico
LS <- qnorm(0.975,mt,sqrt(Ct))
LI <- qnorm(0.025,mt,sqrt(Ct))

plot(yt)
lines(dropFirst(ft), col = 2)
lines(LS, col = 4, lty = 2)
lines(LI, col = 4, lty = 2)
