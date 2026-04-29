library(readr)
library(tseries)
library(forecast)
library(urca)               
library(lmtest)
library(nortest)
library(VGAM)
library(zoo)
library(mass) 

### Lectura y visualización inicial ####

#filepath = "C:/Users/yetli/Documents/UDLAP/PR 2026/Eco II/Time-Series/Final/New/Consumo.csv"
filepath = "D:/Gu/GEN/Eco/Final/New/Consumo.csv"
df =  read.csv(filepath, fileEncoding = "latin1")
head(df)
df
consumoTS <- ts(df$Consumo, 
                start = c(1993, 1),
                end = c(2025, 12),
                frequency = 12)
consumoTS
?ts
plot(consumoTS, type = "l", 
     main = "Indicador mensual del consumo privado", 
     sub = "Base 2018. Serie desestacionalizada",xlab = "Fecha", 
     ylab = "Índice de volumen físico 2018 = 100", lwd = 2)

ggseasonplot(consumoTS) #no es estacional :)

### Estacionariedad ####

adf.test(consumoTS) # p-value = 0.2316
ndiffs(consumoTS) # d=1
nsdiffs(consumoTS) # d=0

# Original
consumoTS_D1 = diff(consumoTS, differences=1)
adf.test(consumoTS_D1) # p-value = 0.01
var(consumoTS_D1) # 1.492427

# Box Cox
lambdaTS = BoxCox.lambda(consumoTS) # -0.9999242
consumoTS_BC = BoxCox(consumoTS, lambdaTS)
adf.test(consumoTS_BC) #p-value = 0.5287

consumoTS_BC_D1 = diff(consumoTS_BC, differences=1)
adf.test(consumoTS_BC_D1) #p-value = 0.01
var(consumoTS_BC_D1) #3.395698e-08

### Valores atípicos ####
time_values <- time(consumoTS)

outliers_OG <- tsoutliers(consumoTS)

plot(consumoTS, type = "l", 
     main = "Serie Original con Valores Atípicos", 
     xlab = "Tiempo", ylab = "Sin Diferenciar")

# Use time_values to map positions → real x-axis coordinates
points(time_values[outliers_OG$index], 
       consumoTS[outliers_OG$index], 
       col = "red", pch = 19, cex = 1.5)


### FAC y FACP ####

par(mfrow = c(1, 1))
# Original
acf(consumoTS, main = "Autocorrelación Simple - Serie Original", 
    xlab="Retrasos")   
pacf(consumoTS, main = "Autocorrelación Parcial - Serie Original", 
     xlab="Retrasos", ylab="PACF")

# Original D1
acf(consumoTS_D1, main = "Autocorrelación Simple - Serie Diferenciada d=1", 
    xlab="Retrasos") # MA = 1
pacf(consumoTS_D1, main = "Autocorrelación Parcial - Serie Diferenciada d=1", 
     xlab="Retrasos", ylab="PACF") # AR = 0 OR 1

# Box Cox D1
acf(consumoTS_BC_D1, main = "Autocorrelación Simple - Serie Transformada Box Cox d=1", 
    xlab="Retrasos") # MA = 1
pacf(consumoTS_BC_D1, main = "Autocorrelación Parcial - Serie Transformada Box Cox d=1", 
     xlab="Retrasos", ylab="PACF") # AR = 0




### Modelos propuestos ####

# Modelo 1: Auto arima
m1_AUTO = auto.arima(consumoTS)

# Modelo 2: Original d=1 ARIMA(0,1,1)
m2_D1_011 

# Modelo 3: Original d=1 ARIMA(1,1,1)
m3_D1_111

# Modelo 4: Box Cox d=1 ARIMA(0,1,1)
m4_BC_011
