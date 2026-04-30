### Libraries ####
library(readr)
library(tseries)
library(forecast)
library(urca)               
library(lmtest)
library(nortest)
library(zoo)
library(tsoutliers)
library(VGAM)

### Lectura y visualización inicial ####

filepath = "C:/Users/yetli/Documents/UDLAP/PR 2026/Eco II/Time-Series/Final/New/Consumo.csv"
#filepath = "D:/Gu/GEN/Eco/Final/New/Consumo.csv"
df =  read.csv(filepath, fileEncoding = "latin1")
consumoTS <- ts(df$Consumo, 
                start = c(1993, 1),
                end = c(2025, 12),
                frequency = 12)

plot(consumoTS, type = "l", 
     main = "Indicador mensual del consumo privado", 
     sub = "Base 2018. Serie desestacionalizada",xlab = "Fecha", 
     ylab = "Índice de volumen físico 2018 = 100", lwd = 2)

ggseasonplot(consumoTS,
             main="Gráfico estacional: Indicador mensual del consumo privado",
             xlab="Mes",
             ylab="Índice de volumen físico 2018 = 100") #no es estacional :)

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

par(mfrow = c(1, 1))

axis_tiempo <- time(consumoTS)

outliers_OG <- tsoutliers(consumoTS)
plot(consumoTS, type = "l", 
     main = "Serie Original con Valores Atípicos", 
     xlab = "Tiempo", ylab = "Sin Diferenciar")
points(axis_tiempo[outliers_OG$index], 
       consumoTS[outliers_OG$index], 
       col = "red", pch = 19, cex = 1.5)


outliers_BC <- tsoutliers(consumoTS_BC)
plot(consumoTS_BC, type = "l", 
     main = "Serie Transformada con Valores Atípicos", 
     xlab = "Tiempo", ylab = "Sin Diferenciar")
points(axis_tiempo[outliers_BC$index], 
       consumoTS_BC[outliers_BC$index], 
       col = "red", pch = 19, cex = 1.5)

### FAC y FACP ####

par(mfrow = c(2, 1))
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

# Modelo 1: Auto arima ARIMA(2,1,2)
m1_AUTO = auto.arima(consumoTS)

# Modelo 2: Original d=1 ARIMA(0,1,1)
m2_D1_011 = Arima(consumoTS, order=c(0,1,1))

# Modelo 3: Original d=1 ARIMA(1,1,1)
m3_D1_111 = arima(consumoTS, order=c(1,1,1))

# Modelo 4: Box Cox d=1 ARIMA(0,1,1)
m4_BC_011 = arima(consumoTS_BC, order=c(0,1,1))

### Verificación de supuestos ####
# Valores atípicos e influyentes ####
tso()
?cook.distance
# Errores media cero #####
residuos_m1 = residuals(m1_AUTO)
residuos_m2 = residuals(m2_D1_011)
residuos_m3 = residuals(m3_D1_111)
residuos_m4 = residuals(m4_BC_011)

?t.test

t.test(residuos_m1) # p.value: 0.9936443
t.test(residuos_m2) # p.value: 0.03929234
t.test(residuos_m3) # p.value: 0.0363425
t.test(residuos_m4) # p.value: 0.009259041

# Errores homocedásticos ####

fitted_m1 = fitted(m1_AUTO)
fitted_m2 = fitted(m2_D1_011)
fitted_m3 = fitted(m3_D1_111)
fitted_m4 = fitted(m4_BC_011)

bptest(residuos_m1^2 ~ fitted_m1) # p-value = 0.5062
bptest(residuos_m2^2 ~ fitted_m2) # p-value = 0.5023
bptest(residuos_m3^2 ~ fitted_m3) # p-value = 0.5002
bptest(residuos_m4^2 ~ fitted_m4) # p-value = 0.593

# Errores normales ####

ad.test(residuos_m1)
ad.test(residuos_m2)
ad.test(residuos_m3)
ad.test(residuos_m4)

shapiro.test(residuos_m1)
shapiro.test(residuos_m2)
shapiro.test(residuos_m3)
shapiro.test(residuos_m4)

lillie.test(residuos_m1)
lillie.test(residuos_m2)
lillie.test(residuos_m3)
lillie.test(residuos_m4)

qqnorm(residuos_m1,
       main = "QQ Plot - Residuos Modelo 1: Auto ARIMA")
qqline(residuos_m1, col = "red")
qqnorm(residuos_m2,
       main = "QQ Plot - Residuos Modelo 2: ARIMA(0,1,1)")
qqline(residuos_m2, col = "red")
qqnorm(residuos_m3,
       main = "QQ Plot - Residuos Modelo 3: ARIMA(1,1,1)")
qqline(residuos_m3, col = "red")
qqnorm(residuos_m4,
       main = "QQ Plot - Residuos Modelo 4: Box Cox ARIMA(0,1,1)")
qqline(residuos_m4, col = "red")

hist(residuos_m1, main = "Histograma - Residuos Modelo 1: Auto ARIMA", 
         xlab = "Residuos", breaks = 20, col = "lightblue")
hist(residuos_m2, main = "Histograma - Residuos Modelo 2: ARIMA(0,1,1)",
         xlab = "Residuos", breaks = 20, col = "lightblue")
hist(residuos_m3, main = "Histograma - Residuos Modelo 3: ARIMA(1,1,1)",
         xlab = "Residuos", breaks = 20, col = "lightblue")
hist(residuos_m4, main = "Histograma - Residuos Modelo 4: Box Cox ARIMA(0,1,1)",
         xlab = "Residuos", breaks = 20, col = "lightblue")

          
          
          
# Errores independientes ####

Box.test(residuos_m1, lag=12, type = "Box-Pierce", fitdf = "4")
Box.test(residuos_m2, lag=12, type = "Box-Pierce", fitdf = "1")
Box.test(residuos_m3, lag=12, type = "Box-Pierce", fitdf = "2")
Box.test(residuos_m4, lag=12, type = "Box-Pierce", fitdf = "1")

acf(residuos_m1, main = "ACF - Residuos Modelo 1: Auto ARIMA", xlab = "Retrasos")
acf(residuos_m2, main = "ACF - Residuos Modelo 2: ARIMA(0,1,1)", xlab = "Retrasos")
acf(residuos_m3, main = "ACF - Residuos Modelo 3: ARIMA(1,1,1)", xlab = "Retrasos")
acf(residuos_m4, main = "ACF - Residuos Modelo 4: Box Cox ARIMA(0,1,1)", xlab = "Retrasos")


