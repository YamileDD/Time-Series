# Relaciones no lineales

Aunque la **autocovarianza** y la **autocorrelación** miden relaciones lineales en series de tiempo, también existen herramientas para medir relaciones **no lineales**. Aquí tienes algunas:

---

### **1. Coeficiente de Autocorrelación No Lineal**
Algunas extensiones del coeficiente de autocorrelación tradicional buscan capturar relaciones más complejas. Por ejemplo:
- **Autocorrelación basada en energía**: Relaciona la distancia energética entre puntos en lugar de su linealidad.
- **Autocorrelación de Pearson transformada**: Usa transformaciones no lineales en los datos antes de calcular el coeficiente.

---

### **2. Métodos de Dependencia No Lineal**
1. **Autocovarianza Generalizada**:
   - Considera relaciones no lineales usando funciones transformadas, como \(g(X_t)\) y \(h(X_{t-k})\), en lugar de \(X_t\) directamente.
   - Fórmula:
     $$
     \text{Cov}(g(X_t), h(X_{t-k})).
     $$
     Ejemplo: Usar \(g(x) = x^2\) y \(h(x) = \log(x)\) para medir relaciones cuadráticas o logarítmicas.

2. **Coeficiente de Correlación de Spearman o Kendall**:
   - Miden relaciones **monotónicas** en lugar de estrictamente lineales.
   - Útil si los datos tienen dependencias crecientes o decrecientes no lineales.

---

### **3. Pruebas de Independencia No Lineal**
1. **Prueba de BDS (Brock-Dechert-Scheinkman)**:
   - Detecta dependencias no lineales en los residuos de una serie de tiempo.
   - Utiliza el concepto de dimensiones embebidas y distribuciones conjuntas para determinar si los datos siguen un comportamiento caótico o independiente.

2. **Prueba de Heterocedasticidad ARCH/GARCH**:
   - Evalúa si las dependencias no lineales están presentes en la varianza de la serie.

---

### **4. Modelos No Lineales**
1. **Modelos de Series Temporales No Lineales**:
   - **TAR (Threshold Autoregressive)**:
     - Permite diferentes comportamientos según si la serie cruza un umbral.
   - **STAR (Smooth Transition Autoregressive)**:
     - Los cambios en el comportamiento dependen suavemente de un parámetro.
   - **GARCH (Generalized Autoregressive Conditional Heteroskedasticity)**:
     - Modela la varianza condicional, útil para datos financieros.
   
2. **Modelos Basados en Redes Neuronales**:
   - Capturan relaciones no lineales complejas utilizando aprendizaje profundo.
   - Ejemplo: LSTM (Long Short-Term Memory), que es especialmente útil para dependencias largas y no lineales.

---

### **5. Métodos Basados en Información**
1. **Entropía Mutua (Mutual Information)**:
   - Mide la dependencia estadística entre \(X_t\) y \(X_{t-k}\), no necesariamente lineal.
   - Se basa en la teoría de la información y evalúa cuánta información sobre \(X_t\) está contenida en \(X_{t-k}\).

2. **Transferencia de Entropía (Transfer Entropy)**:
   - Evalúa el flujo de información direccional entre dos series de tiempo.

---

### **6. Análisis Espectral No Lineal**
- Herramientas como la **transformada wavelet** pueden descomponer una serie de tiempo en diferentes frecuencias y detectar patrones no lineales en diferentes escalas temporales.

---

### **Resumen**
Si sospechas que hay **dependencias no lineales** en tus datos, podrías usar herramientas como:
- **Entropía mutua** para relaciones no lineales generales.
- **Modelos TAR/STAR** si hay cambios en el comportamiento según ciertos umbrales.
- **Pruebas ARCH/GARCH** si hay relaciones no lineales en la varianza.
- **Pruebas de BDS** si buscas dependencias caóticas.