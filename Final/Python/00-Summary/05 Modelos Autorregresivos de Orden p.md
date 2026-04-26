# **Modelos Autorregresivos de Orden \(p\) (AR(p))**

Los modelos **autorregresivos** son una de las clases más comunes de modelos para series de tiempo. Se basan en la idea de que los valores futuros de una serie de tiempo pueden ser explicados como una combinación lineal de sus valores pasados.

---

### **Definición del Modelo AR(p)**

Un modelo AR(p) describe una serie de tiempo \(X_t\) como una combinación lineal de sus \(p\) valores pasados y un término de error:

$$
X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \dots + \phi_p X_{t-p} + \epsilon_t,
$$

donde:
- \(X_t\): Valor de la serie en el tiempo \(t\).
- \(\phi_1, \phi_2, \dots, \phi_p\): Coeficientes autorregresivos que determinan el peso de los valores pasados.
- \(p\): Orden del modelo, es decir, el número de valores pasados utilizados.
- \(\epsilon_t\): Error aleatorio o término de ruido blanco, $\epsilon_t \sim N(0, \sigma^2)$.

---

### **Propiedades del Modelo AR(p)**

1. **Dependencia Temporal**:
   - Cada valor de la serie depende linealmente de \(p\) valores anteriores.
   - Los coeficientes \(\phi_i\) reflejan la importancia de cada lag.

2. **Estacionariedad**:
   - Para que el modelo AR(p) sea estacionario, las raíces del polinomio característico:
     $$
     1 - \phi_1 z - \phi_2 z^2 - \dots - \phi_p z^p = 0
     $$
     deben estar fuera del círculo unitario (\(|z| > 1\)).

3. **ACF y PACF**:
   - La función de autocorrelación (ACF) tiende a decaer gradualmente (exponencial o sinusoidal).
   - La función de autocorrelación parcial (PACF) tiene \(p\) términos significativos y se corta después del lag \(p\).

---

### **Casos Particulares**

1. **Modelo AR(1)**:
   $$
   X_t = \phi_1 X_{t-1} + \epsilon_t.
   $$
   - Dependencia solo del valor inmediatamente anterior.
   - La ACF decae exponencialmente.

2. **Modelo AR(2)**:
   $$
   X_t = \phi_1 X_{t-1} + \phi_2 X_{t-2} + \epsilon_t.
   $$
   - Dependencia de los dos valores anteriores.
   - Puede generar comportamientos más complejos, como oscilaciones (si \(\phi_1\) y \(\phi_2\) están cerca de ciertos valores).

---

### **Estimación de Parámetros**
Para estimar los coeficientes \(\phi_1, \phi_2, \dots, \phi_p\), se utilizan métodos como:
1. **Máxima Verosimilitud**:
   - Maximiza la probabilidad de observar los datos dados los parámetros del modelo.
2. **Método de Yule-Walker**:
   - Resuelve un sistema de ecuaciones basado en las autocorrelaciones de la serie.
3. **Mínimos Cuadrados Ordinarios (OLS)**:
   - Ajusta el modelo minimizando el error cuadrático.

---

### **Diagnóstico del Modelo**
1. **Validación de Residuos**:
   - Los residuos deben ser ruido blanco (independientes y con media cero).
   - Se verifica con la prueba de Ljung-Box o analizando la ACF de los residuos.

2. **AIC y BIC**:
   - Se utilizan criterios de información como el AIC (Criterio de Información de Akaike) y BIC (Criterio de Información Bayesiano) para seleccionar el orden \(p\).

---

### **Ejemplo Práctico**
Supongamos una serie de tiempo \(X_t\) que sigue un modelo AR(2):
$$
X_t = 0.5 X_{t-1} - 0.2 X_{t-2} + \epsilon_t.
$$

1. Si observas los valores \(X_{t-1} = 2.0\) y \(X_{t-2} = 1.5\), y el error es \(\epsilon_t = 0.3\), entonces:
   $$
   X_t = 0.5(2.0) - 0.2(1.5) + 0.3 = 1.0 - 0.3 + 0.3 = 1.0.
   $$

2. Su ACF decaerá gradualmente con posibles oscilaciones.

---

### **Aplicaciones del Modelo AR(p)**
1. **Finanzas**: Modelar precios de activos financieros o tasas de interés.
2. **Economía**: Predecir indicadores económicos como PIB, inflación o ventas.
3. **Meteorología**: Predecir temperaturas diarias o niveles de precipitación.