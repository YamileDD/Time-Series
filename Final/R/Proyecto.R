install.packages("tseries")
install.packages("forecast")
install.packages("urca")
install.packages("lmtest")
install.packages("nortest")
install.packages("VGAM")

library(readr)
library(tseries)
library(forecast)
library(urca)               
library(lmtest)
library(nortest)
library(VGAM)

datos = read.csv("opcion3.csv", fileEncoding = "latin1")
datos_indice = datos$Índice

# Gráfico de la serie original
plot(datos_indice, type = "l", main = "Serie Temporal Original", 
     xlab = "Tiempo", ylab = "Valor del Índice", col = "hotpink", lwd = 2)

# Identificación automática de modelo ARIMA
auto.arima(datos_indice)
length(datos_indice) # Número de observaciones


# VERIFICACIÓN DE ESTACIONARIEDAD

# Prueba de raíz unitaria (ADF)
adf.test(datos_indice)
# Conclusión: No es estacionario

# Análisis gráfico de correlaciones serie original
acf(datos_indice, main = "Autocorrelación - Serie Original", xlab="Retrasos")   
pacf(datos_indice, main = "Autocorrelación Parcial - Serie Original", xlab="Retrasos", ylab="PACF") 

# Aplicación de primera diferencia para lograr estacionariedad
dif1_bc = diff(datos_indice)

# Verificación de estacionariedad en serie diferenciada
adf.test(dif1_bc)
# Conclusión: Con una diferencia se vuelve estacionario

# Análisis de correlaciones serie diferenciada
acf(dif1_bc, main = "Autocorrelación - Serie Diferenciada", xlab="Retrasos")  
pacf(dif1_bc, main = "Autocorrelación Parcial - Serie Diferenciada", xlab="Retrasos", ylab="PACF")  
plot(dif1_bc, type = "l", main = "Serie Temporal Diferenciada", 
     xlab = "Tiempo", ylab = "Primera Diferencia")


# MODELO INICIAL COMPLEJO: ARIMA(5,1,3)

ARIMA513 = arima(datos_indice, order = c(5,1,3))

# Prueba de media cero en residuos
t.test(resid(ARIMA513), conf.level = 0.95)
# P-val > 0.05, no rechazamos H0, los residuos tienen media 0

# Estandarización de residuos para análisis
res_ARIMA513 = resid(ARIMA513)
res_ARIMA513_std = res_ARIMA513/sqrt(ARIMA513$sigma2)
plot(res_ARIMA513_std, type = "l", main = "Residuales Estandarizados", 
     xlab = "Tiempo", ylab = "Residuales Estandarizados")
abline(h = c(-2, 2), col = "red", lty = 2)

# Detección de valores atípicos
tsoutliers(dif1_bc)

# Gráfico de residuales del modelo
plot(resid(ARIMA513), main = "Residuales del Modelo ARIMA(5,1,3)", 
     xlab = "Tiempo", ylab = "Residuales")

# Prueba de homocedasticidad - Breusch-Pagan
mod513 = lm((resid(ARIMA513))^2 ~ fitted(ARIMA513))
bptest(mod513)

# Histograma de residuales
hist(resid(ARIMA513), breaks = 50, main = "Histograma de Residuales", 
     xlab = "Residuales", ylab = "Frecuencia")

# Prueba de normalidad - Gráfico Q-Q y Anderson-Darling
qqnorm(resid(ARIMA513), main = "Gráfico Q-Q de Residuales", xlab="Cuantiles Teóricos", ylab="Cuantiles Muestrales")
qqline(resid(ARIMA513))
ad.test(resid(ARIMA513))
# P-val pequeño, rechazamos H0, los residuos no siguen normalidad

# Gráfico con límites de control para residuales
desv_d1 = sqrt(var(resid(ARIMA513)))
LI2D_d1 = -2*desv_d1  # Límite inferior 2 desviaciones
LS2D_d1 = 2*desv_d1   # Límite superior 2 desviaciones
plot(resid(ARIMA513), main = "Residuos con Límites de Control", 
     ylim = c(-10, 8), xlab = "Tiempo", ylab = "Residuales")
abline(h = LI2D_d1, col = "hotpink")
abline(h = LS2D_d1, col = "hotpink")

# Prueba de independencia - Box-Pierce
Box.test(resid(ARIMA513), lag = 1, type = "Box-Pierce")

# Intervalos de confianza para parámetros
confint(ARIMA513, level = 0.95) 

# Matriz de varianza-covarianza
ARIMA513$var.coef
cov2cor(ARIMA513$var.coef)  # Matriz de correlaciones

# Verificación de estabilidad - Raíces de polinomios
ar_roots = polyroot(c(1, -ARIMA513$coef[grep("ar", names(ARIMA513$coef))]))
Mod(ar_roots)  # Módulo de raíces AR
ma_roots = polyroot(c(1, ARIMA513$coef[grep("ma", names(ARIMA513$coef))]))
Mod(ma_roots)  # Módulo de raíces MA


# MODELO 1: ARIMA(1, 1, 0) - MODELO SELECCIONADO

ARIMA110 = arima(datos_indice, order = c(1,1,0))
ARIMA110

# Verificación de supuestos para ARIMA(1,1,0)
t.test(resid(ARIMA110), conf.level = 0.95)
# P-val grande, no rechazamos H0, los residuos tienen media 0

# Homocedasticidad
mod110 = lm((resid(ARIMA110))^2 ~ fitted(ARIMA110))
bptest(mod110)
# P-val grande, hay homocedasticidad

plot(resid(ARIMA110), main = "Residuales del Modelo ARIMA(1,1,0)", 
     xlab = "Tiempo", ylab = "Residuales")
abline(h = 0, lty = 2, col = "red", lwd = 2)

# Normalidad en los errores
hist(resid(ARIMA110), breaks = 50, col = "hotpink", 
     main = "Histograma de Residuales - ARIMA(1,1,0)", 
     xlab = "Residuales", ylab = "Frecuencia")
qqnorm(resid(ARIMA110), main = "Gráfico Q-Q - ARIMA(1,1,0)", xlab="Cuantiles Teóricos", ylab="Cuantiles Muestrales")
qqline(resid(ARIMA110))
ad.test(resid(ARIMA110))
# P-val muy pequeño, rechazamos H0, residuos no siguen normalidad

# Errores independientes
Box.test(resid(ARIMA110), lag = 1, type = "Box-Pierce")
acf(residuals(ARIMA110), main = "Autocorrelación de Residuales-ARIMA(1,1,0)", xlab="Retrasos")
pacf(residuals(ARIMA110), main = "Autocorrelación Parcial de Residuales-ARIMA(1,1,0)", xlab="Retrasos", ylab="PACF")

# Parsimonia
confint(ARIMA110, level = 0.99) 

# Estabilidad del modelo
ar_roots_110 = polyroot(c(1, -ARIMA110$coef[grep("ar", names(ARIMA110$coef))]))
Mod(ar_roots_110)

# Valores atípicos
tsoutliers(dif1_bc)
plot(dif1_bc, type = "l", main = "Serie Diferenciada con Valores Atípicos", 
     xlab = "Tiempo", ylab = "Primera Diferencia")
points(tsoutliers(dif1_bc)$index, 
       dif1_bc[tsoutliers(dif1_bc)$index], 
       col="red", pch=19)


# MODELO 2: ARIMA(0, 1, 1) CON DRIFT

ARIMA011 = Arima(datos_indice, order = c(0,1,1), include.drift = TRUE)
ARIMA011
# Verificación de supuestos para ARIMA(0,1,1)
t.test(resid(ARIMA011), conf.level = 0.95)
# P-val grande, no rechazamos H0, los residuos tienen media 0

# Homocedasticidad
mod011 = lm((resid(ARIMA011))^2 ~ fitted(ARIMA011))
bptest(mod011)
# P-val grande, hay homocedasticidad

plot(resid(ARIMA011), main = "Residuales del Modelo ARIMA(0,1,1)", 
     xlab = "Tiempo", ylab = "Residuales")
abline(h = 0, lty = 2, col = "red", lwd = 2)

# Normalidad en los errores
hist(resid(ARIMA011), breaks = 50, col = "hotpink", 
     main = "Histograma de Residuales - ARIMA(0,1,1)", 
     xlab = "Residuales", ylab = "Frecuencia")
qqnorm(resid(ARIMA011), pch = 1, main = "Gráfico Q-Q - ARIMA(0,1,1)", xlab="Cuantiles Teóricos", ylab="Cuantiles Muestrales")
qqline(resid(ARIMA011), col = "red", lwd = 2)
ad.test(resid(ARIMA011))
# P-val muy pequeño, rechazamos H0, residuos no siguen normalidad

# Errores independientes
Box.test(resid(ARIMA011), lag = 1, type = "Box-Pierce")
acf(residuals(ARIMA011), main = "Autocorrelación de Residuales - ARIMA(0,1,1)", xlab="Retrasos")
pacf(residuals(ARIMA011), main = "Autocorrelación Parcial de Residuales - ARIMA(0,1,1)", xlab="Retrasos", ylab="PACF")

# Parsimonia
confint(ARIMA011, level = 0.99) 

# Estabilidad del modelo
ma_roots_011 = polyroot(c(1, ARIMA011$coef[grep("ma", names(ARIMA011$coef))]))
Mod(ma_roots_011)


# MODELO 3: ARIMA(1, 1, 1)

ARIMA111 = arima(datos_indice, order = c(1,1,1))
ARIMA111

# Verificación de supuestos para ARIMA(1,1,1)
t.test(resid(ARIMA111), conf.level = 0.95)
# P-val grande, no rechazamos H0, los residuos tienen media 0

# Homocedasticidad
mod111 = lm((resid(ARIMA111))^2 ~ fitted(ARIMA111))
bptest(mod111)
# P-val grande, hay homocedasticidad

plot(resid(ARIMA111), main = "Residuales del Modelo ARIMA(1,1,1)", 
     xlab = "Tiempo", ylab = "Residuales")
abline(h = 0, lty = 2, col = "red", lwd = 2)

# Normalidad en los errores
hist(resid(ARIMA111), breaks = 50, col = "hotpink", 
     main = "Histograma de Residuales - ARIMA(1,1,1)", 
     xlab = "Residuales", ylab = "Frecuencia")
qqnorm(resid(ARIMA111), pch = 1, main = "Gráfico Q-Q - ARIMA(1,1,1)", xlab="Cuantiles Teóricos", ylab="Cuantiles Muestrales")
qqline(resid(ARIMA111), col = "red", lwd = 2)
ad.test(resid(ARIMA111))
# P-val muy pequeño, rechazamos H0, residuos no siguen normalidad

# Errores independientes
Box.test(resid(ARIMA111), lag = 1, type = "Box-Pierce")
acf(residuals(ARIMA111), main = "Autocorrelación de Residuales - ARIMA(1,1,1)", xlab="Retrasos")
pacf(residuals(ARIMA111), main = "Autocorrelación Parcial de Residuales - ARIMA(1,1,1)", xlab="Retrasos", ylab="PACF")

# Parsimonia
confint(ARIMA111, level = 0.99) 

# Estabilidad del modelo
ar_roots_111 = polyroot(c(1, -ARIMA111$coef[grep("ar", names(ARIMA111$coef))]))
Mod(ar_roots_111)
ma_roots_111 = polyroot(c(1, ARIMA111$coef[grep("ma", names(ARIMA111$coef))]))
Mod(ma_roots_111)

cat("COMPARACIÓN DE MODELOS ARIMA:\n")
cat("ARIMA(1,1,0) - AIC:", AIC(ARIMA110), "BIC:", BIC(ARIMA110), "\n")
cat("ARIMA(0,1,1) con drift - AIC:", AIC(ARIMA011), "BIC:", BIC(ARIMA011), "\n")
cat("ARIMA(1,1,1) - AIC:", AIC(ARIMA111), "BIC:", BIC(ARIMA111), "\n\n")

# CÁLCULO DE PRONÓSTICOS PARA MODELO ARIMA(1,1,0) SELECCIONADO
  
modelo_final = ARIMA110
phi = as.numeric(modelo_final$coef["ar1"])
sigma = sqrt(modelo_final$sigma2)

cat("PARÁMETROS DEL MODELO ARIMA(1,1,0) SELECCIONADO\n")
cat("φ =", phi, "\n")
cat("σ_a =", sigma, "\n\n")

# FUNCIÓN PARA CÁLCULO DE PRONÓSTICOS
calc_pred = function(phi, x_t, x_t1, h_max) {
  pred = numeric(h_max)
  
  for(h in 1:h_max) {
    if(h == 1) {
      pred[h] = (1 + phi) * x_t - phi * x_t1
    } else if(h == 2) {
      pred[h] = (1 + phi) * pred[1] - phi * x_t
    } else {
      pred[h] = (1 + phi) * pred[h-1] - phi * pred[h-2]
    }
  }
  return(pred)
}

# FUNCIÓN PARA CÁLCULO DE COEFICIENTES ψ
calc_psi = function(phi, h_max) {
  psi = numeric(h_max + 1)
  psi[1] = 1
  psi[2] = 1 + phi
  
  for(j in 3:(h_max + 1)) {
    psi[j] = (1 + phi) * psi[j-1] - phi * psi[j-2]
  }
  return(psi)
}

# FUNCIÓN PARA INTERVALOS DE PREDICCIÓN
calc_int = function(psi, pred, h_max, sigma) {
  z = qnorm(0.975)  # Valor crítico para 95% de confianza
  int = matrix(0, nrow = h_max, ncol = 3)
  
  for(h in 1:h_max) {
    var_error = sigma^2 * sum(psi[1:h]^2)
    error_estandar = sqrt(var_error)
    LI = pred[h] - z * error_estandar
    LS = pred[h] + z * error_estandar
    int[h,] = c(LI, pred[h], LS)
  }
  return(int)
}

# PRONÓSTICOS
ultimos = tail(datos_indice, 2)
x_t = ultimos[1]    # Penúltimo valor de la serie
x_t1 = ultimos[2]   # Último valor de la serie
cat("PRONÓSTICOS SERIE ORIGINAL (10 periodos):\n")
pred_orig = calc_pred(phi, x_t, x_t1, 10)
for(i in 1:10) {
  cat("Pronóstico", i, ":", pred_orig[i], "\n")
}

cat("\nINTERVALOS DE PREDICCIÓN SERIE ORIGINAL (95% confianza):\n")
psi_coef = calc_psi(phi, 10)
int_orig = calc_int(psi_coef, pred_orig, 10, sigma)
for(i in 1:10) {
  cat("Periodo", i, ": [", round(int_orig[i,1], 4), ", ", round(int_orig[i,3], 4), "]\n")
}
# ACTUALIZACIÓN DE PRONÓSTICOS CON NUEVA OBSERVACIÓN
nuevo_valor = 117.4546204
innov = nuevo_valor - ((1 + phi) * x_t - phi * x_t1)
nuevo_x_t = nuevo_valor
nuevo_x_t1 = x_t

cat("\nPRONÓSTICOS ACTUALIZADOS (9 periodos):\n")
pred_act = calc_pred(phi, nuevo_x_t, nuevo_x_t1, 9)
for(i in 1:9) {
  cat("Pronóstico", i, ":", pred_act[i], "\n")
}

cat("\nINTERVALOS DE PREDICCIÓN ACTUALIZADOS (95% confianza):\n")
int_act = calc_int(psi_coef, pred_act, 9, sigma)
for(i in 1:9) {
  cat("Periodo", i, ": [", round(int_act[i,1], 4), ", ", round(int_act[i,3], 4), "]\n")
}
# VERIFICACIÓN FINAL CON FUNCIÓN FORECAST
fc = forecast(modelo_final, h = 10)
comp = data.frame(
  Manual = round(pred_orig, 6),
  Forecast = round(as.numeric(fc$mean), 6),
  Diferencia = round(pred_orig - as.numeric(fc$mean), 8)
)

# Resultados de la comparación
comp

