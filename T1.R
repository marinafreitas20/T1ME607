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
library(kableExtra)

#Função para plotar FAC, FACP e qqplot dos resíduos
residuos_plots <- function (fit,lag ){
  acf(fit$residuals, lag.max = lag)
  pacf(fit$residuals, lag.max = lag)
  qqnorm(fit$residuals)
  qqline(fit$residuals, col = "red")
}


#Função para mostar o p-valor para vários testes de normalidade dos resíduos
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


# Limpando e orgizando o banco de dados -----------------------------------

#Leitura dos dados
dados <- read.table('Chicago-65ate74.txt')
#Renomeando as colunas dos dados
colnames(dados)<-c("Data","Mortes")
#transformando a variável data em formato de data
dados$Data <- as.Date(dados$Data, format ='%Y-%m-%d')

#Quantos NAs temos nos dados
sum(is.na(dados$Mortes))
#Qual o índice de cada um dos NAs
which(is.na(dados$Mortes))

#Substituindo NAs pela média dos vizinhos
dados$Mortes[54] = (dados$Mortes[53]+dados$Mortes[55])/2
dados$Mortes[3118] = (dados$Mortes[3117]+dados$Mortes[3121])/2
dados$Mortes[3119] = (dados$Mortes[3118]+dados$Mortes[3121])/2 #Já usamos a média anterior
dados$Mortes[3120] = (dados$Mortes[3119]+dados$Mortes[3121])/2 #Já usamos a média anterior

#Encontrando qual dia da semana a data pertence
dados$diasem <- weekdays(dados$Data) 

#Numerando a qual semana cada data pertence
dados$sem <- c(1,1,1,rep(c(2:731),  each = 7),732)

#Novo data frame com os óbitos agrupados por semana
dados_sem <- dados %>% group_by(sem) %>% summarise(mortes = sum(Mortes))

#Deixamos na base a data onde inicia a semana
dados_sem$Data <- c(dados$Data[1],dados$Data[dados$diasem == 'domingo'])


# Questão 1, item A -------------------------------------------------------

############################# Plotando o gráfico da série para ver como ela
##### se comporta e quais são suas carterísticas empíricas

ggplot(dados_sem, aes(x = sem, y = mortes))+
  geom_line()+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+
  theme_classic()


############################# Plotando o gráfico da FAC e FACP da série 
##### para ver como elas se comportam
par(mfrow=c(1,1))

acf(as.ts(dados_sem$mortes), lag=700)
pacf(dados_sem$mortes,  lag=700)


# Questão 1, item B -------------------------------------------------------

############################# Removendo as 13 últimas observações (semanas) 
##### da série

nrow(dados_sem)-13

dados_serie <- dados_sem[1:(nrow(dados_sem)-13),]

##Transformando os dados sem as ultimas 13 obs em serie temporal
#Decidimos tirar também a primeira semana (semana 1) e a última (semana 732), já que os valores não correspondem a semanas completas
serie = ts(dados_serie$mortes[2:718], start=c(1987,1), frequency = 52)
#### Os dados para testar o modelo escolhido será então dados_serie$mortes[719:731]

####################################################################################
################################## Harmônico #######################################
####################################################################################

tam <- length(serie)        # tamanho da serie
t.1 <- seq(1:tam)           # tendencia
n.p <- tam/52                # numero de periodos
nu <- 1/52

####################################### Harmônico 1 ################################
#Com tendência quadrática e sazonalidade com um termo de seno e um termo de cosseno

#tendência quadrática
t.2 <- t.1^2

preditoras_harmonico <- data.frame(t.1, t.2, nu = rep(nu,717))

#Ajuste do modelo
ajuste.h1 <- lm(serie ~ t.1 + t.2 + sin(2*pi*nu*t.1) + cos(2*pi*nu*t.1), data = preditoras_harmonico) 
#Coeficientes do modelo
summary(ajuste.h1)

#Dados para plotagem da série original com a série ajustada
df_serie_1 <- data.frame(dados_serie$sem[2:718],serie,ajuste.h1$fitted.values )
colnames(df_serie_1) <- c('Semana','observado', 'ajustado')

#Plotando série + ajuste
ggplot(df_serie_1, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red", size =1.25)+
  theme_classic()+
  ylab("Óbitos")

#Teste dos residuos
residuos_plots(ajuste.h1, lag = 350)
residuos_testes(ajuste.h1, lag = 350)
Box.test(ajuste.h1$residuals, type = "Ljung-Box")

#Valores de AIC e BIC
AIC(ajuste.h1)
BIC(ajuste.h1)


BL2 = function(y,p,dmax=20){
  v = c()
  for (i in (p+1):dmax){v[i] = Box.test(y, lag=i-p, type="Ljung-Box")$p.value}
  plot((p+1):dmax,v[(p+1):dmax],pch=20,ylim=c(0,1),xlab="Defasagem",ylab="Valor-p")
  abline(h=0.05,lty=2,col="blue")
}



########################################################################################
############################ DUMMY #####################################################
########################################################################################



d <- seasonaldummy(serie)

for(i in 1:13){
  d[52*i,] <- rep(-1,51)  }        

preditoras <- data.frame(t.1,d)

ajuste.d <- lm(serie ~ t.1 + (S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+S21+
                                S22+S23+S24+S25+S26+S27+S28+S29+S30+S31+S32+S33+S34+S35+S36+S37+S38+S39+S40+
                                S41+S42+S43+S44+S45+S46+S47+S48+S49+S50+S51 ), data = preditoras)


summary(ajuste.d)

df_serie_dummy <- data.frame(dados_serie$sem[2:718],serie[2:718],ajuste.d$fitted.values[2:718] )
colnames(df_serie_dummy) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =0.7)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d)
residuos_testes(ajuste.d)
Box.test(ajuste.d$residuals, type = "Ljung-Box")
ljung

AIC(ajuste.d)
BIC(ajuste.d)

####################################################################################
############################## PREDICT #############################################
####################################################################################

######################### Dados do teste ###########################################

teste$mortes <- dados_sem$mortes[719:731]
teste$semana <- week(dados_sem$Data[719:731])


############################## Predict do modelo Dummy ############################

t.teste <- c(719:731)
dummy_teste <- d[c(40:52),]

df.teste_dummy <- data.frame(t.1 = t.teste,dummy_teste)

forecast_dummy <- forecast(object = ajuste.d, newdata = df.teste_dummy)

forecast_dummy

df_forecast_dummy <- data.frame(c(719:731),dados_sem$mortes[719:731],forecast_dummy$mean )
colnames(df_forecast_dummy) <- c('Semana','observado', 'ajustado')

ggplot(df_forecast_dummy, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =0.7)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Óbitos" )+ 
  theme_classic()

#Erros
dados_sem$mortes[719:731]-forecast_dummy$mean

#Soma dos erros quadrados
sum((dados_sem$mortes[719:731]-forecast_dummy$mean)^2)

############################## Predict do modelo Harmônico ############################

t.teste <- c(719:731)
t2.teste <- t.teste^2

df.teste_harmonico <- data.frame(t.1 = t.teste,t.2 = t2.teste,nu = rep(nu,13))

forecast_harmonico <- forecast(object = ajuste.h1, newdata = df.teste_harmonico)

forecast_harmonico

df_forecast_harmonico <- data.frame(c(719:731),dados_sem$mortes[719:731],forecast_harmonico$mean)
colnames(df_forecast_harmonico) <- c('Semana','observado', 'ajustado')

ggplot(df_forecast_harmonico, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =0.7)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Óbitos" )+ 
  theme_classic()

#Erros
dados_sem$mortes[719:731] - forecast_harmonico$mean

#Soma dos erros quadrados
sum((dados_sem$mortes[719:731] - forecast_harmonico$mean)^2)

#******************************************************************************************
#******************************TENTATIVAS DE OUTROS MODELOS *******************************
#******************************************************************************************

########################################################################################
############################ DUMMY  com tendência quadrática ###########################
########################################################################################

t.2 <- t.1^2

preditoras_quadra <- data.frame(t.1,t.2,d)

ajuste.d.quadratico <- lm(serie ~ t.1 + t.2 + (S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+S21+
                                                 S22+S23+S24+S25+S26+S27+S28+S29+S30+S31+S32+S33+S34+S35+S36+S37+S38+S39+S40+
                                                 S41+S42+S43+S44+S45+S46+S47+S48+S49+S50+S51 ), data = preditoras_quadra)

summary(ajuste.d.quadratico )

df_serie_dummy_quadratico <- data.frame(dados_serie$sem[2:718],serie[2:718],ajuste.d.quadratico$fitted.values[2:718])
colnames(df_serie_dummy_quadratico) <- c('Semana','observado', 'ajustado')

ggplot(df_serie_dummy_quadratico, aes(x=Semana)) + 
  geom_line(aes(y = observado)) + 
  geom_line(aes(y = ajustado), color="red",  size =0.8)+
  scale_x_continuous(breaks = c(0:14)*52) +
  xlab("Semana") + ylab("Número de óbitos" )+ 
  theme_classic()

#Teste dos residuos
residuos_plots(ajuste.d.quadratico )
residuos_testes(ajuste.d.quadratico )
Box.test(ajuste.d.quadratico $residuals, type = "Ljung-Box")

AIC(ajuste.d.quadratico)
BIC(ajuste.d.quadratico)

################################################################################
############## Harmônico com 2 senos e cossenos ################################
################################################################################

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

###################################################################################
########################## Séries transformadas #####################################
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
############################ DUMMY com sqrt ############################################
########################################################################################


ajuste.d.sqrt <- lm(sqrt(serie) ~ t.1 + (S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+S21+
                                           S22+S23+S24+S25+S26+S27+S28+S29+S30+S31+S32+S33+S34+S35+S36+S37+S38+S39+S40+
                                           S41+S42+S43+S44+S45+S46+S47+S48+S49+S50+S51 ), data = preditoras)

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


ajuste.d.ln <- lm(log(serie) ~ t.1 + (S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+S21+
                                        S22+S23+S24+S25+S26+S27+S28+S29+S30+S31+S32+S33+S34+S35+S36+S37+S38+S39+S40+
                                        S41+S42+S43+S44+S45+S46+S47+S48+S49+S50+S51 ), data = preditoras)

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


ajuste.d.Box <- lm(BoxCox(serie, lambda = "auto") ~ t.1 + (S1+S2+S3+S4+S5+S6+S7+S8+S9+S10+S11+S12+S13+S14+S15+S16+S17+S18+S19+S20+S21+
                                                             S22+S23+S24+S25+S26+S27+S28+S29+S30+S31+S32+S33+S34+S35+S36+S37+S38+S39+S40+
                                                             S41+S42+S43+S44+S45+S46+S47+S48+S49+S50+S51 ), data = preditoras)

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




