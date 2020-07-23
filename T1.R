setwd("G:/Meu Drive/Unicamp/ME607/Trabalho 1/T1ME607")

library(tidyverse)
library(lubridate)
library(ggplot2)
library(lmtest)
library(nortest)
library(forecast)
library(astsa)
library(lattice)
library(timeSeries)

residuos_plots <- function (fit){
  acf(fit$residuals, lag.max = 100)
  pacf(fit$residuals, lag.max = 100)
  qqnorm(fit$residuals)
  qqline(fit$residuals, col = "red")
}



residuos_testes <- function (fit) { 
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


acf(as.ts(dados_serie$mortes),main = 'Mortes Chicago', lag=700  )
pacf(dados_serie$mortes, main = 'Mortes Chicago', lag=700)

serie = ts(dados_serie$mortes[2:718], start=c(1987,1), frequency = 52) ##Transformando os dados sem as ultimas 13 obs em serie temporal

aditivo <-  decompose(serie, type = 'additive')

plot(aditivo)

multiplicativo <-  decompose(serie, type = 'multiplicative')

plot(multiplicativo)


###########################################################
################### ARMA ##################################
###########################################################

############### ARMA nos dados originais ##################

auto.arima(dados_serie$mortes)

fit_arma <- arima(serie, order = c(1,0,1))
summary(fit_arma)

df_serie_ARMA <- data.frame(dados_serie$sem[2:718],serie,fitted(fit_arma) )
colnames(df_serie_ARMA) <- c('Semana','observado', 'ajustado')

#Plotando série + ajuste
ggplot(df_serie_ARMA, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) +
  theme_classic()


#testes dos resíduos
residuos_plots(fit_arma)
residuos_testes(fit_arma)
Box.test(fit_arma$residuals, type = "Ljung-Box")

############### ARMA com transformação sqrt #################

fit_arma_sqrt <- arima(sqrt(serie), order = c(1,0,1))
summary(fit_arma_sqrt)

df_serie_ARMA_sqrt <- data.frame(dados_serie$sem[2:718],sqrt(serie),fitted(fit_arma_sqrt) )
colnames(df_serie_ARMA_sqrt) <- c('Semana','observado', 'ajustado')

#Plotando série + ajuste
ggplot(df_serie_ARMA_sqrt, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) +
  theme_classic()

#testes dos resíduos
residuos_plots(fit_arma_sqrt)
residuos_testes(fit_arma_sqrt)
Box.test(fit_arma_sqrt$residuals, type = "Ljung-Box")

############### ARMA com transformação logaritmica #################

fit_arma_ln <- arima(log(serie), order = c(1,0,1))
summary(fit_arma_ln)

df_serie_ARMA_ln <- data.frame(dados_serie$sem[2:718],log(serie),fitted(fit_arma_ln) )
colnames(df_serie_ARMA_ln) <- c('Semana','observado', 'ajustado')

#Plotando série + ajuste
ggplot(df_serie_ARMA_ln, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) +
  theme_classic()

#testes dos resíduos
residuos_plots(fit_arma_ln)
residuos_testes(fit_arma_ln)
Box.test(fit_arma_ln$residuals, type = "Ljung-Box")


###########################################################
################### Harmônico #############################
###########################################################

tam <- length(serie)        # tamanho da serie
t.1 <- seq(1:tam)           # tendencia
n.p <- tam/52                # numero de periodos
nu <- 1/52

############## Harmônico 1 ################################

ajuste.h1 <- lm(serie ~ t.1 + t.1^2 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1)) 
summary(ajuste.h1)

df_serie_1 <- data.frame(dados_serie$sem[2:718],serie,ajuste.h1$fitted.values )
colnames(df_serie_1) <- c('Semana','observado', 'ajustado')

#Plotando série + ajuste
ggplot(df_serie_1, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.h1)
residuos_testes(ajuste.h1)
Box.test(ajuste.h1$residuals, type = "Ljung-Box")


###################################################################################
########################## Série transformada #####################################
###################################################################################

############## Harmônico transformada sqrt ################################

ajuste.hs <- lm(sqrt(serie) ~ t.1 + t.1^2+ sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) )

df_serie_hs <- data.frame(dados_serie$sem[2:718],sqrt(serie),ajuste.hs$fitted.values )
colnames(df_serie_hs) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_hs, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.hs)
residuos_testes(ajuste.hs)
Box.test(ajuste.hs$residuals, type = "Ljung-Box")

############## Harmônico transformada ln ################################

ajuste.hln <- lm(log(serie) ~ t.1 + t.1^2+ sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) )

df_serie_hln <- data.frame(dados_serie$sem[2:718],log(serie),ajuste.hln$fitted.values )
colnames(df_serie_hln) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_hln, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.hln)
residuos_testes(ajuste.hln)
Box.test(ajuste.hln$residuals, type = "Ljung-Box")

############## Harmônico transformada Box Cox ################################

ajuste.hBox <- lm(BoxCox(serie, lambda = "auto") ~ t.1 + t.1^2+ sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) )

df_serie_hBox <- data.frame(dados_sem$sem[2:718],BoxCox(serie, lambda = "auto"),ajuste.hBox$fitted.values )
colnames(df_serie_hBox) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_hBox, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.hBox)
residuos_testes(ajuste.hBox)
Box.test(ajuste.hBox$residuals, type = "Ljung-Box")

################################################################################
############## Harmônico com 2 senos e cossenos ################################

ajuste.h2 <- lm(serie ~ t.1 + t.1^2 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) + sin(2*pi*2*nu*t.1) + cos(2*pi*2*nu*t.1))

df_serie_h2 <- data.frame(dados_serie$sem[2:718],serie,ajuste.h2$fitted.values )
colnames(df_serie_h2) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_h2, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.h2)
residuos_testes(ajuste.h2)
Box.test(ajuste.h2$residuals, type = "Ljung-Box")

############## Harmônico 2 com transformação sqrt ################################

ajuste.h2_sqrt <- lm(sqrt(serie) ~ t.1 + t.1^2 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) + sin(2*pi*2*nu*t.1) + cos(2*pi*2*nu*t.1))

df_serie_h2_sqrt <- data.frame(dados_serie$sem[2:718],sqrt(serie),ajuste.h2_sqrt$fitted.values )
colnames(df_serie_h2_sqrt) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_h2_sqrt, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.h2_sqrt)
residuos_testes(ajuste.h2_sqrt)
Box.test(ajuste.h2_sqrt$residuals, type = "Ljung-Box")

############## Harmônico 2 com transformação logaritmica ################################

ajuste.h2_ln <- lm(log(serie) ~ t.1 + t.1^2 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) + sin(2*pi*2*nu*t.1) + cos(2*pi*2*nu*t.1))

df_serie_h2_ln <- data.frame(dados_serie$sem[2:718],log(serie),ajuste.h2_ln$fitted.values )
colnames(df_serie_h2_ln) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_h2_ln, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.h2_ln)
residuos_testes(ajuste.h2_ln)
Box.test(ajuste.h2_ln$residuals, type = "Ljung-Box")

############## Harmônico 2 com transformação de Box Cox ################################

ajuste.h2_Box <- lm(BoxCox(serie, lambda = "auto") ~ t.1 + t.1^2 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1) + sin(2*pi*2*nu*t.1) + cos(2*pi*2*nu*t.1))

df_serie_h2_Box <- data.frame(dados_serie$sem[2:718],BoxCox(serie, lambda = "auto"),ajuste.h2_Box$fitted.values )
colnames(df_serie_h2_Box) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_h2_Box, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", linetype="twodash", size =2) 

#Teste dos residuos
residuos_plots(ajuste.h2_Box)
residuos_testes(ajuste.h2_Box)
Box.test(ajuste.h2_Box$residuals, type = "Ljung-Box")

########################################################################################
############################ DUMMY #####################################################
########################################################################################



d <- seasonaldummy(serie)

for(i in 1:13){
  d[52*i,] <- rep(-1,51)  }        


ajuste.d <- lm(serie ~ t.1 + (d[,1]+d[,2]+d[,3]+d[,4]+d[,5]+d[,6]+d[,7]+d[,8]+d[,9]+d[,10]+
                                     d[,11]+d[,12]+d[,13]+d[,14]+d[,15]+d[,16]+d[,17]+d[,18]+d[,19]+d[,20]+
                                     d[,21]+d[,22]+d[,23]+d[,24]+d[,25]+d[,26]+d[,27]+d[,28]+d[,29]+d[,30]+
                                     d[,31]+d[,32]+d[,33]+d[,34]+d[,35]+d[,36]+d[,37]+d[,38]+d[,39]+d[,40]+
                                     d[,41]+d[,42]+d[,43]+d[,44]+d[,45]+d[,46]+d[,47]+d[,48]+d[,49]+d[,50]+
                                     d[,51] ))

summary(ajuste.d)

df_serie_dummy <- data.frame(dados_serie$sem[2:718],serie,ajuste.d$fitted.values )
colnames(df_serie_dummy) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =1.25)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d)
residuos_testes(ajuste.d)
Box.test(ajuste.d$residuals, type = "Ljung-Box")

AIC(ajuste.d)
BIC(ajuste.d)

########################################################################################
############################ DUMMY  QUADRATICO #########################################
########################################################################################

t.2 <- t.1^2

ajuste.d.quadratico <- lm(serie ~ t.1 + t.2 + (d[,1]+d[,2]+d[,3]+d[,4]+d[,5]+d[,6]+d[,7]+d[,8]+d[,9]+d[,10]+
                                d[,11]+d[,12]+d[,13]+d[,14]+d[,15]+d[,16]+d[,17]+d[,18]+d[,19]+d[,20]+
                                d[,21]+d[,22]+d[,23]+d[,24]+d[,25]+d[,26]+d[,27]+d[,28]+d[,29]+d[,30]+
                                d[,31]+d[,32]+d[,33]+d[,34]+d[,35]+d[,36]+d[,37]+d[,38]+d[,39]+d[,40]+
                                d[,41]+d[,42]+d[,43]+d[,44]+d[,45]+d[,46]+d[,47]+d[,48]+d[,49]+d[,50]+
                                d[,51] ))

summary(ajuste.d.quadratico )

df_serie_dummy_quadratico <- data.frame(dados_serie$sem[2:718],serie,ajuste.d.quadratico$fitted.values )
colnames(df_serie_dummy_quadratico) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy_quadratico, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =1.25)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d.quadratico )
residuos_testes(ajuste.d.quadratico )
Box.test(ajuste.d.quadratico $residuals, type = "Ljung-Box")

AIC(ajuste.d.quadratico)
BIC(ajuste.d.quadratico)

########################################################################################
############################ DUMMY com sqrt ############################################
########################################################################################


ajuste.d.sqrt <- lm(sqrt(serie) ~ t.1 + (d[,1]+d[,2]+d[,3]+d[,4]+d[,5]+d[,6]+d[,7]+d[,8]+d[,9]+d[,10]+
                                d[,11]+d[,12]+d[,13]+d[,14]+d[,15]+d[,16]+d[,17]+d[,18]+d[,19]+d[,20]+
                                d[,21]+d[,22]+d[,23]+d[,24]+d[,25]+d[,26]+d[,27]+d[,28]+d[,29]+d[,30]+
                                d[,31]+d[,32]+d[,33]+d[,34]+d[,35]+d[,36]+d[,37]+d[,38]+d[,39]+d[,40]+
                                d[,41]+d[,42]+d[,43]+d[,44]+d[,45]+d[,46]+d[,47]+d[,48]+d[,49]+d[,50]+
                                d[,51] ))

summary(ajuste.d.sqrt)

df_serie_dummy_sqrt <- data.frame(dados_serie$sem[2:718],sqrt(serie),ajuste.d.sqrt$fitted.values )
colnames(df_serie_dummy_sqrt) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy_sqrt, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =1.25)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d.sqrt)
residuos_testes(ajuste.d.sqrt)
Box.test(ajuste.d.sqrt$residuals, type = "Ljung-Box")

AIC(ajuste.d.sqrt)
BIC(ajuste.d.sqrt)

########################################################################################
############################ DUMMY com transformação logaritmica #######################
########################################################################################


ajuste.d.ln <- lm(log(serie) ~ t.1 + (d[,1]+d[,2]+d[,3]+d[,4]+d[,5]+d[,6]+d[,7]+d[,8]+d[,9]+d[,10]+
                                           d[,11]+d[,12]+d[,13]+d[,14]+d[,15]+d[,16]+d[,17]+d[,18]+d[,19]+d[,20]+
                                           d[,21]+d[,22]+d[,23]+d[,24]+d[,25]+d[,26]+d[,27]+d[,28]+d[,29]+d[,30]+
                                           d[,31]+d[,32]+d[,33]+d[,34]+d[,35]+d[,36]+d[,37]+d[,38]+d[,39]+d[,40]+
                                           d[,41]+d[,42]+d[,43]+d[,44]+d[,45]+d[,46]+d[,47]+d[,48]+d[,49]+d[,50]+
                                           d[,51] ))

summary(ajuste.d.ln)

df_serie_dummy_ln <- data.frame(dados_serie$sem[2:718],log(serie),ajuste.d.ln$fitted.values )
colnames(df_serie_dummy_ln) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy_ln, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =1.25)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d.ln)
residuos_testes(ajuste.d.ln)
Box.test(ajuste.d.ln$residuals, type = "Ljung-Box")

AIC(ajuste.d.ln)
BIC(ajuste.d.ln)

########################################################################################
############################ DUMMY com transformação de Box Cox ########################
########################################################################################


ajuste.d.Box <- lm(BoxCox(serie, lambda = "auto") ~ t.1 + (d[,1]+d[,2]+d[,3]+d[,4]+d[,5]+d[,6]+d[,7]+d[,8]+d[,9]+d[,10]+
                                        d[,11]+d[,12]+d[,13]+d[,14]+d[,15]+d[,16]+d[,17]+d[,18]+d[,19]+d[,20]+
                                        d[,21]+d[,22]+d[,23]+d[,24]+d[,25]+d[,26]+d[,27]+d[,28]+d[,29]+d[,30]+
                                        d[,31]+d[,32]+d[,33]+d[,34]+d[,35]+d[,36]+d[,37]+d[,38]+d[,39]+d[,40]+
                                        d[,41]+d[,42]+d[,43]+d[,44]+d[,45]+d[,46]+d[,47]+d[,48]+d[,49]+d[,50]+
                                        d[,51] ))

summary(ajuste.d.Box)

df_serie_dummy_Box <- data.frame(dados_serie$sem[2:718],BoxCox(serie, lambda = "auto"),ajuste.d.Box$fitted.values )
colnames(df_serie_dummy_Box) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy_Box, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =1.25)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d.Box)
residuos_testes(ajuste.d.Box)
Box.test(ajuste.d.Box$residuals, type = "Ljung-Box")

AIC(ajuste.d.Box)
BIC(ajuste.d.Box)


####################################################################################
############################## PREDICT #############################################
####################################################################################

teste$mortes <- dados_sem$mortes[718:731]
teste$semana <- week(dados_sem$Data[718:731])

dummy_teste <- d[c(39:52),]
t.teste <- c(718:731)

df.teste <- data.frame(t.teste,dummy_teste)

predict.lm(object = ajuste.d, newdata = df.teste)

