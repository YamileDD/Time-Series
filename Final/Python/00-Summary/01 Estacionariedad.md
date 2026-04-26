# **Estacionariedad**

### **¿Qué es una serie de tiempo?**
Es un conjunto de datos ordenados en el tiempo, donde el objetivo principal es analizar su comportamiento y predecir valores futuros con base en patrones históricos.

### **Definición de Estacionariedad**
Una serie de tiempo es **estacionaria** si sus propiedades estadísticas, como la media, la varianza y la covarianza, son constantes a lo largo del tiempo. Esto implica que el proceso subyacente que genera la serie no cambia.

#### **Propiedades clave de una serie estacionaria:**
1. **Media constante**: No varía con el tiempo.
   - Ejemplo: $E[X_t] = \mu$, donde $\mu$ es constante.
2. **Varianza constante**: No depende del tiempo.
   - Ejemplo: $\text{Var}(X_t) = \sigma^2$.
3. **Covarianza constante**: La relación entre los valores de la serie depende solo de la distancia entre ellos (\(k\)), no del tiempo específico.
   - Ejemplo: $\text{Cov}(X_t, X_{t-k}) = \gamma_k$, donde $\gamma_k$ depende solo de $k$.

---

### **Importancia de la Estacionariedad**
- **Modelos más simples**: Muchos métodos econométricos, como ARMA y ARIMA, requieren estacionariedad para realizar inferencias válidas.
- **Predicciones confiables**: Sin estacionariedad, las relaciones estadísticas entre los datos no se mantienen a lo largo del tiempo.
- **Evitar regresiones espurias**: Si se analiza una serie no estacionaria sin transformarla, los resultados pueden ser engañosos.

---

### **Tipos de Estacionariedad**
1. **Estacionariedad estricta**:
   - Todas las distribuciones conjuntas (de cualquier orden) son invariantes en el tiempo.
   - Es una definición teórica y difícil de cumplir en la práctica.
2. **Estacionariedad débil (o en media y varianza)**:
   - Solo requiere que la media, varianza y covarianza sean constantes.
   - Es más práctica y la más usada en econometría.

---

### **Cómo detectar estacionariedad**
1. **Análisis gráfico**:
   - Observa si la serie tiene una tendencia o si la varianza parece cambiar con el tiempo.
2. **Pruebas estadísticas**:
   - **Prueba de Dickey-Fuller Aumentada (ADF)**:
     - Hipótesis nula (\(H_0\)): La serie tiene una raíz unitaria (no es estacionaria).
     - Hipótesis alternativa (\(H_1\)): La serie es estacionaria.
   - **Prueba KPSS (Kwiatkowski-Phillips-Schmidt-Shin)**:
     - \(H_0\): La serie es estacionaria.
     - \(H_1\): La serie no es estacionaria.
3. **Correlograma (ACF y PACF)**:
   - La ACF (Función de autocorrelación) decae rápidamente en una serie estacionaria.

---

### **Cómo transformar una serie para que sea estacionaria**
1. **Diferenciación**:
   - Substraer el valor de un periodo anterior: $Y_t' = Y_t - Y_{t-1}$.
   - Si una diferencia no es suficiente, se puede aplicar varias veces (diferenciación de orden \(d\)).
2. **Logaritmos**:
   - Útil para estabilizar la varianza.
3. **Descomposición de tendencia**:
   - Separar la tendencia y analizar los residuos.
4. **Filtros**:
   - Usar filtros como Hodrick-Prescott para eliminar ciclos o tendencias.

---

### **Ejemplo práctico:**
Supongamos que tenemos una serie de tiempo del PIB trimestral de un país. 

1. Observamos un gráfico y notamos una **tendencia ascendente**. Esto sugiere que **no es estacionaria**.
2. Aplicamos una **diferenciación de primer orden**:
   - $Y_t' = Y_t - Y_{t-1}$.
3. Verificamos con la prueba ADF y encontramos que después de la diferenciación, \(H_0\) es rechazada.
4. Concluimos que la serie diferenciada es **estacionaria** y podemos usarla para modelar.


