# Cargar librerías y datos 
paquetes <- c("readr", "zoo", "tseries", "lubridate", "forecast", "lmtest")

instalar <- paquetes[!(paquetes %in% installed.packages()[,"Package"])]

if(length(instalar) > 0){
  install.packages(instalar)
}

lapply(paquetes, library, character.only = TRUE)

serie <- read_csv("C:/Users/HUAWEI/Downloads/serie.csv",
                  show_col_types = FALSE)

# Transformación de fechas y serie de tiempo 
serie$Periodos <- as.yearmon(serie$Periodos, format = "%Y/%m")

inicio_anio <- as.numeric(format(min(serie$Periodos), "%Y"))
inicio_mes  <- as.numeric(format(min(serie$Periodos), "%m"))

serie_ts <- ts(
  serie$BIE,
  start = c(inicio_anio, inicio_mes),
  frequency = 12
)

graphics.off()
par(mar = c(4,4,2,1))

plot(serie_ts,
     main = "Serie original",
     col = "blue")

# Prueba de estacionariedad (ADF) -----------------------------------

cat("PRUEBA ADF - SERIE ORIGINAL\n")
resultado_adf <- adf.test(serie_ts)
print(resultado_adf)

# Transformación estabilizadora de varianza --------------------------

lambda_bc <- BoxCox.lambda(serie_ts)

cat("TRANSFORMACIÓN BOX-COX\n")
cat("Lambda estimado:", lambda_bc, "\n")

serie_bc <- BoxCox(serie_ts, lambda = lambda_bc)

plot(serie_bc,
     main = paste("Serie transformada Box-Cox, lambda =", round(lambda_bc, 4)),
     col = "purple")

# Uso del operador diferencia -----------------------------------

cat("ADF - SERIE TRANSFORMADA\n")
resultado_adf_bc <- adf.test(serie_bc)
print(resultado_adf_bc)

# ACF y PACF --------------------------------------------------------------

datos_est <- serie_bc

par(mfrow = c(1,2))

acf(datos_est,
    main = "ACF de la serie transformada",
    lag.max = 24)

pacf(datos_est,
     main = "PACF de la serie transformada",
     lag.max = 24)

par(mfrow = c(1,1))

# Modelos propuestos  -----------------------------------------------------
# ARIMA (2,0,0)
# ARIMA (0,0,2)
# ARIMA (1,0,2)

modelo_200 <- Arima(datos_est,
                    order = c(2, 0, 0),
                    include.mean = TRUE)

modelo_002 <- Arima(datos_est,
                    order = c(0, 0, 2),
                    include.mean = TRUE)

modelo_120 <- Arima(datos_est,
                    order = c(1, 0, 2),
                    include.mean = TRUE)

# Resumen de cada modelo 

summary(modelo_200)
summary(modelo_002)
summary(modelo_120)

# Estimadores calculados en R

modelo_200$coef
modelo_002$coef
modelo_120$coef

# Comparación de modelos 

AIC(modelo_200, modelo_002, modelo_120)
BIC(modelo_200, modelo_002, modelo_120)

# Verificacion de supuestos  ---------------------------------------------

res_200 <- na.omit(as.numeric(residuals(modelo_200)))
fit_200 <- na.omit(as.numeric(fitted(modelo_200)))

res_002 <- na.omit(as.numeric(residuals(modelo_002)))
fit_002 <- na.omit(as.numeric(fitted(modelo_002)))

res_120 <- na.omit(as.numeric(residuals(modelo_120)))
fit_120 <- na.omit(as.numeric(fitted(modelo_120)))

# Supuesto 1: Errores con media cero --------------------------------------

cat("\nSUPUESTO 1: MEDIA CERO\n")

cat("\nModelo ARIMA(2,0,0)\n")
t.test(res_200, mu = 0)

cat("\nModelo ARIMA(0,0,2)\n")
t.test(res_002, mu = 0)

cat("\nModelo ARIMA(1,0,2)\n")
t.test(res_120, mu = 0)

# Supuesto 2: Varianza constante ------------------------------------------

cat("\nSUPUESTO 2: VARIANZA CONSTANTE\n")

cat("\nModelo ARIMA(2,0,0)\n")
plot(fit_200, res_200,
     main = "Residuos vs Ajustados ARIMA(2,0,0)",
     xlab = "Ajustados",
     ylab = "Residuos")
abline(h = 0, col = "red")
lmtest::bptest(lm(res_200^2 ~ fit_200))

cat("\nModelo ARIMA(0,0,2)\n")
plot(fit_002, res_002,
     main = "Residuos vs Ajustados ARIMA(0,0,2)",
     xlab = "Ajustados",
     ylab = "Residuos")
abline(h = 0, col = "red")
lmtest::bptest(lm(res_002^2 ~ fit_002))

cat("\nModelo ARIMA(1,0,2)\n")
plot(fit_120, res_120,
     main = "Residuos vs Ajustados ARIMA(1,0,2)",
     xlab = "Ajustados",
     ylab = "Residuos")
abline(h = 0, col = "red")
lmtest::bptest(lm(res_120^2 ~ fit_120))

# Supuesto 3: Normalidad de los errores -----------------------------------

cat("\nSUPUESTO 3: NORMALIDAD\n")

cat("\nModelo ARIMA(2,0,0)\n")
shapiro.test(res_200)
hist(res_200, main = "Histograma ARIMA(2,0,0)", col = "gray")
qqnorm(res_200); qqline(res_200, col = "red")

cat("\nModelo ARIMA(0,0,2)\n")
shapiro.test(res_002)
hist(res_002, main = "Histograma ARIMA(0,0,2)", col = "gray")
qqnorm(res_002); qqline(res_002, col = "red")

cat("\nModelo ARIMA(1,0,2)\n")
shapiro.test(res_120)
hist(res_120, main = "Histograma ARIMA(1,0,2)", col = "gray")
qqnorm(res_120); qqline(res_120, col = "red")

# Supuesto 4: Independencia -----------------------------------------------

cat("\nSUPUESTO 4: INDEPENDENCIA\n")

cat("\nModelo ARIMA(2,0,0)\n")

Box.test(res_200,
         lag = 12,
         type = "Box-Pierce",
         fitdf = 2)

Box.test(res_200,
         lag = 12,
         type = "Ljung-Box",
         fitdf = 2)

acf(res_200,
    main = "ACF residuos ARIMA(2,0,0)",
    lag.max = 24)

cat("\nModelo ARIMA(0,0,2)\n")

Box.test(res_002,
         lag = 12,
         type = "Box-Pierce",
         fitdf = 2)

Box.test(res_002,
         lag = 12,
         type = "Ljung-Box",
         fitdf = 2)

acf(res_002,
    main = "ACF residuos ARIMA(0,0,2)",
    lag.max = 24)

cat("\nModelo ARIMA(1,0,2)\n")

Box.test(res_120,
         lag = 12,
         type = "Box-Pierce",
         fitdf = 3)

Box.test(res_120,
         lag = 12,
         type = "Ljung-Box",
         fitdf = 3)

acf(res_120,
    main = "ACF residuos ARIMA(1,0,2)",
    lag.max = 24)

# Supuesto 5: Parsimonia --------------------------------------------------

cat("\nSUPUESTO 5: PARSIMONIA\n")

cat("\nModelo ARIMA(2,0,0)\n")
cat("AIC:", AIC(modelo_200), "\n")
cat("BIC:", BIC(modelo_200), "\n")

cat("\nModelo ARIMA(0,0,2)\n")
cat("AIC:", AIC(modelo_002), "\n")
cat("BIC:", BIC(modelo_002), "\n")

cat("\nModelo ARIMA(1,0,2)\n")
cat("AIC:", AIC(modelo_120), "\n")
cat("BIC:", BIC(modelo_120), "\n")

# Supuesto 6: Modelo admisible --------------------------------------------

cat("\nSUPUESTO 6: MODELO ADMISIBLE\n")

summary(modelo_200)
summary(modelo_002)
summary(modelo_120)

# Supuesto 7: Modelo estable ----------------------------------------------

cat("\nSUPUESTO 7: MODELO ESTABLE\n")

ar_200 <- coef(modelo_200)[grep("ar", names(coef(modelo_200)))]
if(length(ar_200) > 0){
  cat("\nRaíces AR ARIMA(2,0,0)\n")
  print(Mod(polyroot(c(1, -ar_200))))
}

ma_002 <- coef(modelo_002)[grep("ma", names(coef(modelo_002)))]
if(length(ma_002) > 0){
  cat("\nRaíces MA ARIMA(0,0,2)\n")
  print(Mod(polyroot(c(1, ma_002))))
}

ar_120 <- coef(modelo_120)[grep("ar", names(coef(modelo_120)))]
ma_120 <- coef(modelo_120)[grep("ma", names(coef(modelo_120)))]

if(length(ar_120) > 0){
  cat("\nRaíces AR ARIMA(1,0,2)\n")
  print(Mod(polyroot(c(1, -ar_120))))
}

if(length(ma_120) > 0){
  cat("\nRaíces MA ARIMA(1,0,2)\n")
  print(Mod(polyroot(c(1, ma_120))))
}

# Supuesto 8: Valores atípicos --------------------------------------------

cat("\nSUPUESTO 8: VALORES ATÍPICOS\n")

res_est_200 <- res_200 / sd(res_200)
res_est_002 <- res_002 / sd(res_002)
res_est_120 <- res_120 / sd(res_120)

cat("\nModelo ARIMA(2,0,0)\n")
print(which(abs(res_est_200) > 3))

cat("\nModelo ARIMA(0,0,2)\n")
print(which(abs(res_est_002) > 3))

cat("\nModelo ARIMA(1,0,2)\n")
print(which(abs(res_est_120) > 3))

# Tabla resumen de supuestos ----------------------------------------------

evaluar_modelo <- function(modelo, res, fit, p, q){
  
  p_media <- t.test(res, mu = 0)$p.value
  p_var   <- lmtest::bptest(lm(res^2 ~ fit))$p.value
  p_norm  <- shapiro.test(res)$p.value
  p_ind   <- Box.test(res, lag = 12, type = "Ljung-Box", fitdf = p + q)$p.value
  
  media_ok <- ifelse(p_media > 0.05, "✔️", "❌")
  var_ok   <- ifelse(p_var > 0.05, "✔️", "❌")
  norm_ok  <- ifelse(p_norm > 0.05, "✔️", "❌")
  ind_ok   <- ifelse(p_ind > 0.05, "✔️", "❌")
  
  aic_val <- round(AIC(modelo), 2)
  
  ar <- coef(modelo)[grep("ar", names(coef(modelo)))]
  estabilidad <- "✔️"
  if(length(ar) > 0){
    raices <- Mod(polyroot(c(1, -ar)))
    if(any(raices <= 1.01)) estabilidad <- "⚠️"
  }
  
  res_est <- res / sd(res)
  n_out <- sum(abs(res_est) > 3)
  outlier <- ifelse(n_out == 0, "✔️",
                    ifelse(n_out <= 2, "⚠️", "❌"))
  
  return(c(media_ok, var_ok, norm_ok, ind_ok, aic_val, estabilidad, outlier))
}

tabla <- data.frame(
  Supuesto = c("Media cero", "Varianza", "Normalidad",
               "Independencia", "AIC", "Estabilidad", "Outliers"),
  
  ARIMA_200 = evaluar_modelo(modelo_200, res_200, fit_200, 2, 0),
  ARIMA_002 = evaluar_modelo(modelo_002, res_002, fit_002, 0, 2),
  ARIMA_120 = evaluar_modelo(modelo_120, res_120, fit_120, 1, 2)
)

tabla

# Pronóstico e intervalos de predicción -----------------------------------

modelo_final <- modelo_200

# Pronóstico a 3 periodos
pron_3 <- forecast(modelo_final, h = 3, level = 95)

# Pronóstico a 10 periodos
pron_10 <- forecast(modelo_final, h = 10, level = 95)

# Tablas en escala original -----------------------------------------------

pron_3_df <- data.frame(
  Periodo = 1:3,
  Pronostico = InvBoxCox(pron_3$mean, lambda_bc),
  LI_95 = InvBoxCox(pron_3$lower[,1], lambda_bc),
  LS_95 = InvBoxCox(pron_3$upper[,1], lambda_bc)
)

pron_10_df <- data.frame(
  Periodo = 1:10,
  Pronostico = InvBoxCox(pron_10$mean, lambda_bc),
  LI_95 = InvBoxCox(pron_10$lower[,1], lambda_bc),
  LS_95 = InvBoxCox(pron_10$upper[,1], lambda_bc)
)

pron_3_df
pron_10_df

plot(pron_3,
     main = "Pronóstico ARIMA(2,0,0) - 3 periodos")

plot(pron_10,
     main = "Pronóstico ARIMA(2,0,0) - 10 periodos")

# Comparación de pronósticos 3 vs 10 --------------------------------------

comparacion_pronosticos <- data.frame(
  Periodo = 1:3,
  Pron_3 = pron_3_df$Pronostico,
  Pron_10 = pron_10_df$Pronostico[1:3],
  LI_3 = pron_3_df$LI_95,
  LS_3 = pron_3_df$LS_95,
  LI_10 = pron_10_df$LI_95[1:3],
  LS_10 = pron_10_df$LS_95[1:3]
)

comparacion_pronosticos

# Actualización del pronóstico --------------------------------------------

# Pronóstico original
pron_original <- forecast(modelo_final, h = 3, level = 95)

# Nuevo dato real después del último dato observado
# IMPORTANTE: sustituye este valor por el dato real cuando lo tengas.
# Mientras tanto, se usa el primer pronóstico como dato observado simulado.
nuevo_dato_original <- InvBoxCox(pron_original$mean[1], lambda_bc)

# Como el modelo está en escala Box-Cox, se transforma el nuevo dato
nuevo_dato_bc <- BoxCox(nuevo_dato_original, lambda_bc)

serie_actualizada <- ts(
  c(as.numeric(datos_est), as.numeric(nuevo_dato_bc)),
  start = start(datos_est),
  frequency = frequency(datos_est)
)

modelo_actualizado <- Arima(serie_actualizada,
                            order = c(2,0,0),
                            include.mean = TRUE)

pron_actualizado <- forecast(modelo_actualizado,
                             h = 3,
                             level = 95)

tabla_actualizacion <- data.frame(
  Periodo = 1:3,
  
  Original = InvBoxCox(pron_original$mean, lambda_bc),
  Actualizado = InvBoxCox(pron_actualizado$mean, lambda_bc),
  
  LI_original = InvBoxCox(pron_original$lower[,1], lambda_bc),
  LS_original = InvBoxCox(pron_original$upper[,1], lambda_bc),
  
  LI_actualizado = InvBoxCox(pron_actualizado$lower[,1], lambda_bc),
  LS_actualizado = InvBoxCox(pron_actualizado$upper[,1], lambda_bc)
)

tabla_actualizacion