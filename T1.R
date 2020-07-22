setwd("G:/Meu Drive/Unicamp/ME607/Trabalho 1/T1ME607")
library(lubridate)
library(timeSeries)
library(forecast)
library(dplyr)

testesmodelo_arima <- function (fit){
  acf(fit$residuals, lag.max = 100)
  pacf(fit$residuals, lag.max = 100)
  qqnorm(fit$residuals)
  qqline(fit$residuals, col = "red")
  coeftest(fit)
}



testesresiduos_arima <- function (fit) { 
  t1 <- ks.test(fit$residuals, "pnorm", mean(fit$residuals), sd(fit$residuals)) # KS
  t2 <- lillie.test(fit$residuals) # Lilliefors
  t3 <- cvm.test(fit$residuals) # Cramér-von Mises
  t4 <- shapiro.test(fit$residuals) # Shapiro-Wilk
  t5 <- sf.test(fit$residuals) # Shapiro-Francia
  t6 <- ad.test(fit$residuals) # Anderson-Darling
  # Tabela de resultados
  testes <- c(t1$method, t2$method, t3$method, t4$method, t5$method,
              t6$method)
  estt <- as.numeric(c(t1$statistic, t2$statistic, t3$statistic,
                       t4$statistic, t5$statistic, t6$statistic))
  valorp <- c(t1$p.value, t2$p.value, t3$p.value, t4$p.value, t5$p.value,
              t6$p.value)
  resultados <- cbind(estt, valorp)
  rownames(resultados) <- testes
  colnames(resultados) <- c("Estatística", "p-valor")
  print(resultados, digits = 4)
  
}

testesresiduos_decompose <- function (fit) { 
  t1 <- ks.test(fit$random, "pnorm", mean(fit$random), sd(fit$random)) # KS
  t2 <- lillie.test(fit$random) # Lilliefors
  t3 <- cvm.test(fit$random) # Cramér-von Mises
  t4 <- shapiro.test(fit$random) # Shapiro-Wilk
  t5 <- sf.test(fit$random) # Shapiro-Francia
  t6 <- ad.test(fit$random) # Anderson-Darling
  # Tabela de resultados
  testes <- c(t1$method, t2$method, t3$method, t4$method, t5$method,
              t6$method)
  estt <- as.numeric(c(t1$statistic, t2$statistic, t3$statistic,
                       t4$statistic, t5$statistic, t6$statistic))
  valorp <- c(t1$p.value, t2$p.value, t3$p.value, t4$p.value, t5$p.value,
              t6$p.value)
  resultados <- cbind(estt, valorp)
  rownames(resultados) <- testes
  colnames(resultados) <- c("Estatística", "p-valor")
  print(resultados, digits = 4)
  
}

testesIsadora <- function (ts,fit){
  
  for (i in 1:20) {
    print(Box.test(fit$residuals, lag = i, type="Ljung-Box")$p.value)
    
  }
  modelagem = ts - fit$residuals
  a = cbind(modelagem, ts)
  xyplot.ts(a, superpose = TRUE)
}

transformacoes <- function(ts){
  transformadas <- list()
  transformadas[[1]] <- sqrt(ts)
  transformadas[[2]] <- log(ts)
  transformadas[[3]] <- BoxCox(ts, lambda = "auto")
  transformadas
}

dados <- read.table('Chicago-65ate74.txt')
colnames(dados)<-c("Data","Mortes")
dados$Data <- as.Date(dados$Data, format ='%Y-%m-%d')


dados$diasem <- weekdays(dados$Data) 

dados$sem <- c(0,0,0,rep(c(1:730),  each = 7),731)

dados_sem <- dados %>% group_by(sem) %>% summarise(mortes = sum(Mortes))

sum(is.na(dados_sem$mortes))

which(is.na(dados_sem$mortes))


dados_sem$mortes[9] = (dados_sem$mortes[8]+dados_sem$mortes[10])/2
dados_sem$mortes[446] = (dados_sem$mortes[445]+dados_sem$mortes[448])/2
dados_sem$mortes[447] = (dados_sem$mortes[446]+dados_sem$mortes[448])/2

dados_sem$Data <- c(dados$Data[1],dados$Data[dados$diasem == 'domingo'])

#ts.plot(dados_sem$mortes, main = "Mortes semanais em Chicago")
#ts.plot(dados_sem$mortes,gpars=list(xlab="Time", ylab= 'mortes', lty=c(1:3)) )


ggplot(dados_sem, aes(x = sem, y = mortes))+
  geom_line()+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+
  theme_classic()
  
nrow(dados_sem)-13

dados_serie <- dados_sem[1:(nrow(dados_sem)-13),]

auto.arima(dados_serie$mortes)

acf(as.ts(dados_serie$mortes),main = 'Mortes Chicago', lag=700  )
pacf(dados_serie$mortes, main = 'Mortes Chicago', lag=700)

serie = ts(dados_serie$mortes, start=c(1987,1), frequency = 52) ##Transformando os dados sem as ultimas 13 obs em serie temporal

serie_transformada <- transformacoes(serie)
####################################################
################ Modelo Aditivo ###################
###################################################

aditivo <-  decompose(serie, type = 'additive')

plot(aditivo)

names(aditivo)
#"x"        "seasonal" "trend"    "random"   "figure"   "type"
aditivo$figure
aditivo$trend
aditivo$seasonal

residual_aditivo <- ts(aditivo$random[27:693], start=c(1987,1), frequency = 52)
acf(residual_aditivo, lag.max = 700)
pacf(residual_aditivo, lag.max = 700)
qqnorm(residual_aditivo)
qqline(residual_aditivo, col = "red")


testesresiduos_decompose(aditivo)
for (i in 1:20) {
  print(Box.test(residual_aditivo, lag = i, type="Ljung-Box")$p.value)
  
}

############ Aditivo raiz quadrada #####################################

aditivo_sqrt <- decompose(serie_transformada[[1]], type = 'additive')

plot(aditivo_sqrt)

residual_aditivo_sqrt <- ts(aditivo_sqrt$random[27:693], start=c(1987,1), frequency = 52)
acf(residual_aditivo_sqrt, lag.max = 700)
pacf(residual_aditivo_sqrt, lag.max = 700)
qqnorm(residual_aditivo_sqrt)
qqline(residual_aditivo_sqrt, col = "red")


testesresiduos_decompose(aditivo_sqrt)
for (i in 1:20) {
  print(Box.test(residual_aditivo_sqrt, lag = i, type="Ljung-Box")$p.value)
  
}

############ Aditivo ln #####################################

aditivo_ln <- decompose(serie_transformada[[2]], type = 'additive')

plot(aditivo_ln)

residual_aditivo_ln <- ts(aditivo_ln$random[27:693], start=c(1987,1), frequency = 52)
acf(residual_aditivo_ln, lag.max = 700)
pacf(residual_aditivo_ln, lag.max = 700)
qqnorm(residual_aditivo_ln)
qqline(residual_aditivo_ln, col = "red")


testesresiduos_decompose(aditivo_ln)
for (i in 1:20) {
  print(Box.test(residual_aditivo_ln, lag = i, type="Ljung-Box")$p.value)
  
}

############ Aditivo Box Cox #####################################

aditivo_Box <- decompose(serie_transformada[[3]], type = 'additive')

plot(aditivo_Box)

residual_aditivo_Box <- ts(aditivo_Box$random[27:693], start=c(1987,1), frequency = 52)
acf(residual_aditivo_Box, lag.max = 700)
pacf(residual_aditivo_Box, lag.max = 700)
qqnorm(residual_aditivo_Box)
qqline(residual_aditivo_Box, col = "red")


testesresiduos_decompose(aditivo_Box)
for (i in 1:20) {
  print(Box.test(residual_aditivo_Box, lag = i, type="Ljung-Box")$p.value)
  
}


####################################################
################ Modelo Multiplicativo #############
###################################################
multiplicativo <- decompose(serie, type = "multiplicative")
plot(multiplicativo)

multiplicativo$figure
multiplicativo$trend
multiplicativo$seasonal

residual_multiplicativo <- ts(multiplicativo$random[27:693], start=c(1987,1), frequency = 52)
acf(residual_multiplicativo, lag.max = 700)
pacf(residual_multiplicativo, lag.max = 700)
qqnorm(residual_multiplicativo)
qqline(residual_multiplicativo, col = "red")


testesresiduos_decompose(multiplicativo)

####################################################
################ Modelo Harmônico #############
###################################################

tam <- length(serie)        # tamanho da serie
t.1 <- seq(1:tam)           # tendencia
n.p <- tam/52                # numero de periodos
nu <- 1/52

ajuste.h1 <- lm(serie ~ t.1 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1)) 
summary(ajuste.h1)
par(mfrow=c(1,1))
ts.plot(serie)
lines(ajuste.h1$fitted.values,lty=4, col='red') # ajuste 

df_serie <- data.frame(dados_sem$sem[1:719],dados_sem$mortes[1:719],ajuste.h1$fitted.values )
colnames(df_serie) <- c('Semana','observado', 'ajustado')

ggplot(df_serie, aes(x=Semana)) + 
  geom_line(aes(y = observado), color = "blue") + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

residual_multiplicativo <- ts(multiplicativo$random[27:693], start=c(1987,1), frequency = 52)
acf(ajuste.h1$residuals, lag.max = 700)
pacf(ajuste.h1$residuals, lag.max = 700)
qqnorm(ajuste.h1$residuals)
qqline(ajuste.h1$residuals, col = "red")

ajuste.h2 <- lm(serie ~ t.1 + t.1^2+ sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) + sin(2*pi*0.025*nu*t.1) + cos(2*pi*0.025*nu*t.1))

df_serie_2 <- data.frame(dados_sem$sem[1:719],dados_sem$mortes[1:719],ajuste.h2$fitted.values )
colnames(df_serie_2) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_2, aes(x=Semana)) + 
  geom_line(aes(y = observado), color = "blue") + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

acf(ajuste.h2$residuals, lag.max = 700)
pacf(ajuste.h2$residuals, lag.max = 700)
qqnorm(ajuste.h2$residuals)
qqline(ajuste.h2$residuals, col = "red")

Box.test(ajuste.h2$residuals, type = "Ljung-Box")
