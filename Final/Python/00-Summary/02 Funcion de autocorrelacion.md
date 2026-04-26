# **Función de autocovarianza**

### **Definición**
La **autocovarianza** mide la relación lineal entre los valores de una serie de tiempo en diferentes momentos, considerando un desfase (o lag) \(k\). Es una medida que cuantifica cuánto se correlacionan dos valores separados por \(k\) periodos.

Matemáticamente, para una serie de tiempo $\{X_t\}$ estacionaria en media y varianza:

$$
\gamma_k = \text{Cov}(X_t, X_{t-k}) = \mathbb{E}[(X_t - \mu)(X_{t-k} - \mu)],
$$
donde:
- $\gamma_k$: Autocovarianza en el lag \(k\).
- $X_t$: Valor de la serie en el tiempo \(t\).
- $\mu$: Media de la serie ($\mathbb{E}[X_t] = \mu$).
- $k$: Desfase o lag.

---

### **Propiedades de la Función de Autocovarianza**
1. **Simetría**: 
   $$
   \gamma_k = \gamma_{-k}.
   $$
   La autocovarianza tiene el mismo valor para un desfase positivo o negativo del mismo tamaño.

2. **Máximo en lag 0**:
   $$
   \gamma_0 = \text{Var}(X_t).
   $$
   En \(k = 0\), la autocovarianza equivale a la varianza de la serie.

3. **Depende solo de \(k\)**:
   - En series estacionarias, la autocovarianza no depende del tiempo \(t\), sino solo del tamaño del desfase \(k\).

---

### **Función de Autocorrelación (ACF)**
Para normalizar la autocovarianza y hacerla más interpretable, se utiliza la **autocorrelación**, que es la autocovarianza dividida por la varianza de la serie:

$$
\rho_k = \frac{\gamma_k}{\gamma_0},
$$
donde $\rho_k$ es la **función de autocorrelación** (ACF) en el lag \(k\).

- $\rho_k \in [-1, 1]$, ya que es un coeficiente de correlación.
- $\rho_k = 1$ significa una relación lineal perfecta.
- $\rho_k = 0$ implica independencia entre los valores en esos momentos.

---

### **Usos de la Función de Autocovarianza y ACF**
1. **Identificar patrones de dependencia temporal**:
   - La forma de la ACF puede revelar si los datos tienen correlaciones persistentes.
2. **Determinar el orden de modelos AR y MA**:
   - En un modelo AR(p), la ACF tiende a decrecer lentamente.
   - En un modelo MA(q), la ACF se corta bruscamente después del lag \(q\).
3. **Evaluar estacionariedad**:
   - En una serie estacionaria, la ACF tiende a decaer rápidamente a medida que aumenta \(k\).

---

### **Gráficos y análisis práctico**
#### **Correlograma**:
- Es un gráfico que muestra la ACF en función del lag \(k\).
- Sirve para visualizar los patrones de dependencia.

#### **Ejemplo:**
Supón que tienes una serie \(X_t = 0.8X_{t-1} + \epsilon_t\), donde $\epsilon_t \sim N(0, \sigma^2)$. 
- La ACF mostrará un decaimiento exponencial, típico de un modelo AR(1).

---



La **diferencia entre correlación y autocorrelación** radica en qué aspectos de los datos están siendo analizados y cómo se mide la relación entre ellos. Aquí te lo explico en detalle:

---

### **Correlación**
- **Definición**: Mide la relación lineal entre dos variables diferentes.
- **Propósito**: Determinar si existe una relación positiva, negativa o nula entre las dos variables.
- **Fórmula (Correlación de Pearson)**:
  $$
  \rho_{XY} = \frac{\text{Cov}(X, Y)}{\sqrt{\text{Var}(X) \cdot \text{Var}(Y)}},
  $$
  donde:
  - $\rho_{XY}$: Correlación entre las variables \(X\) y \(Y\).
  - $\text{Cov}(X, Y)$: Covarianza entre \(X\) y \(Y\).
  - $\text{Var}(X)$ y $\text{Var}(Y)$: Varianzas de \(X\) y \(Y\), respectivamente.

- **Ejemplo**: 
  - Relación entre el ingreso y el gasto en una población.
  - Relación entre la temperatura y el consumo de energía.

---

### **Autocorrelación**
- **Definición**: Mide la relación entre los valores de una misma serie de tiempo en diferentes momentos, con un desfase (\(lag\)).
- **Propósito**: Analizar patrones temporales y dependencias dentro de una serie de tiempo.
- **Fórmula**:
  $$
  \rho_k = \frac{\text{Cov}(X_t, X_{t-k})}{\sqrt{\text{Var}(X_t) \cdot \text{Var}(X_{t-k})}},
  $$
  donde:
  - $\rho_k$: Autocorrelación en el lag \(k\).
  - $\text{Cov}(X_t, X_{t-k})$: Covarianza entre los valores en los tiempos \(t\) y \(t-k\).

- **Ejemplo**:
  - Análisis de la temperatura diaria: ¿Qué tan relacionada está la temperatura de hoy con la de días anteriores?
  - Análisis de ventas mensuales: ¿Qué tan relacionadas están las ventas de un mes con las del mes anterior?

---

### **Diferencias Clave**

| Aspecto              | Correlación                                  | Autocorrelación                              |
|----------------------|----------------------------------------------|---------------------------------------------|
| **Variables**        | Relaciona **dos variables distintas**.       | Relaciona **la misma variable** en diferentes momentos. |
| **Uso**              | Determinar relación lineal entre variables.  | Identificar patrones temporales o dependencias internas. |
| **Dominio**          | Usada en datos transversales o de panel.     | Usada en **series de tiempo**.              |
| **Propósito típico** | Identificar relaciones intervariables.       | Estudiar dependencias temporales.           |

---

### **Ejemplo Comparativo**
- **Correlación**: ¿Cómo se relacionan las ventas de automóviles con el precio del combustible?
- **Autocorrelación**: ¿Cómo se relacionan las ventas de automóviles de este mes con las del mes pasado?