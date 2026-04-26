# **Modelos de Media Móvil de Orden \(q\) (MA(q))**

Los **Modelos de Media Móvil** son otro componente básico en el análisis de series de tiempo. A diferencia de los modelos autorregresivos (AR), los modelos MA utilizan los errores pasados (ruido blanco) para explicar los valores actuales de una serie.

---

### **Definición del Modelo MA(q)**
Un modelo de media móvil de orden \(q\) expresa el valor actual de la serie de tiempo \(X_t\) como una combinación lineal de los errores aleatorios (\(\epsilon_t\)) de los \(q\) periodos anteriores:

$$
X_t = \mu + \epsilon_t + \theta_1 \epsilon_{t-1} + \theta_2 \epsilon_{t-2} + \dots + \theta_q \epsilon_{t-q},
$$

donde:
- \(X_t\): Valor de la serie en el tiempo \(t\).
- \(\mu\): Media de la serie.
- \(\epsilon_t\): Error aleatorio o término de ruido blanco (\(\epsilon_t \sim N(0, \sigma^2)\)).
- \(\theta_1, \theta_2, \dots, \theta_q\): Coeficientes de media móvil que miden la influencia de los errores pasados.
- \(q\): Orden del modelo, es decir, el número de errores pasados utilizados.

---

### **Propiedades del Modelo MA(q)**

1. **No tiene memoria infinita**:
   - En un modelo MA(q), cada valor de \(X_t\) depende únicamente de los errores de los últimos \(q\) periodos.
   - Después de \(q\) periodos, no hay influencia de los errores pasados.

2. **Estacionariedad**:
   - Todo modelo MA(q) es inherentemente estacionario, ya que depende de un número finito de términos de ruido blanco.

3. **ACF y PACF**:
   - La función de autocorrelación (ACF) tiene \(q\) términos significativos y se corta abruptamente después del lag \(q\).
   - La función de autocorrelación parcial (PACF) tiende a decaer gradualmente.

---

### **Interpretación de los Coeficientes**
- Los coeficientes \(\theta_1, \theta_2, \dots, \theta_q\) reflejan cómo los errores pasados (\(\epsilon_{t-1}, \epsilon_{t-2}, \dots\)) afectan el valor actual de la serie.

---

### **Casos Particulares**

1. **Modelo MA(1)**:
   $$
   X_t = \mu + \epsilon_t + \theta_1 \epsilon_{t-1}.
   $$
   - La serie depende únicamente del error actual (\(\epsilon_t\)) y del error inmediato anterior (\(\epsilon_{t-1}\)).

2. **Modelo MA(2)**:
   $$
   X_t = \mu + \epsilon_t + \theta_1 \epsilon_{t-1} + \theta_2 \epsilon_{t-2}.
   $$
   - La serie incorpora la influencia de los dos errores más recientes.

---

### **Estimación de Parámetros**
La estimación de los coeficientes \(\theta_1, \theta_2, \dots, \theta_q\) y la media \(\mu\) se realiza típicamente mediante:
1. **Máxima Verosimilitud**:
   - Encuentra los parámetros que maximizan la probabilidad de los datos observados.
2. **Método de los Momentos**:
   - Relaciona los coeficientes con las autocorrelaciones calculadas.

---

### **Diagnóstico del Modelo**
1. **Validación de Residuos**:
   - Los residuos (\(\epsilon_t\)) deben ser ruido blanco: independientes, con media cero y varianza constante.
2. **ACF**:
   - Verificar si la ACF se corta abruptamente después del lag \(q\), como lo predice el modelo.

---

### **Ejemplo Práctico**
Supongamos que tienes un modelo MA(1):
$$
X_t = 10 + \epsilon_t + 0.5 \epsilon_{t-1}.
$$

1. Si el error actual es \(\epsilon_t = 1.2\) y el error del periodo anterior es \(\epsilon_{t-1} = -0.8\), entonces:
   $$
   X_t = 10 + 1.2 + 0.5(-0.8) = 10 + 1.2 - 0.4 = 10.8.
   $$

2. **ACF**:
   - En este caso, la ACF tendrá un único valor significativo en el lag 1 y será 0 para lags mayores.

---

### **Aplicaciones del Modelo MA(q)**
1. **Suavizado de Series de Tiempo**:
   - Los modelos MA se usan para reducir el ruido en series de tiempo.
2. **Análisis de Residuos**:
   - Los errores en modelos complejos (como ARIMA) suelen seguir un patrón de media móvil.
3. **Modelado en Finanzas**:
   - Utilizados para modelar retornos de activos financieros o volatilidad.

---

### **Relación con el Modelo AR(p)**
| Aspecto             | AR(p)                                      | MA(q)                                      |
|---------------------|--------------------------------------------|--------------------------------------------|
| **Dependencia**     | Basado en valores pasados de la serie.     | Basado en errores pasados.                 |
| **ACF**             | Decae gradualmente.                       | Se corta abruptamente después de \(q\).    |
| **PACF**            | Se corta después de \(p\).                | Decae gradualmente.                        |
| **Estacionariedad** | Puede no ser estacionario.                | Siempre es estacionario.                   |


