import numpy as np

def modulo_raices(**kwargs):
    # Esta función calcula el polinomio característico de un modelo SARIMA y verifica si es estacionario o no.

    # Conocer el coeficiente mas alto
    try:
        max_coef = max([int(i[1:]) for i in kwargs.keys() if i.startswith('p')])
    except: 
        max_coef = 0    

    print("El grado del polinomio es:", max_coef)
    # Inicializar el vector de coeficientes con ceros
    coef = np.zeros(max_coef + 1)
    # Asignar los coeficientes en orden inverso
    for i in range(max_coef, 0, -1):
        coef[max_coef - i] = kwargs.get(f"p{i}", 0)  # x^(i)

    # Agregar el coeficiente de x^0 (que es 1)
    coef[-1] = 1  # x^0

    raices = np.roots(coef)  # Raíces del polinomio
    modulo = np.abs(raices)  # Módulo de las raíces
    estacionario = np.all(modulo > 1)  # Verificar si todas las raíces están fuera del círculo unitario

    print("\nRaíces del polinomio característico:", raices)
    print("\nMódulo de las raíces:", modulo)
    print("\n¿Las raíces están fuera del círculo unitario? ", estacionario)

    return estacionario




def estacionario(modelo):
    "Esta función verifica si el modelo es estacionario o no."
    ars = []
    elevados = []

    # Filtrar solo los parámetros que comienzan con 'ar'
    for i in modelo.params.index:
        if i.startswith('ar'):
            ars.append(-1 * modelo.params[i])   # !!! Se multiplica por -1 
            grado = int(i.split('L')[-1])
            elevados.append(grado)

    # Construir los argumentos para modulo_raices_12.
    # Se asigna cada coeficiente a p{grado} y se dejan en cero los demás.
    params_coef = {}
    for coef, grado in zip(ars, elevados):
        params_coef[f"p{grado}"] = coef

    print("El polinomio de la parte autoregresiva es:", params_coef)
    

    print("")
    # Llamar a la función modulo_raices_12 con los argumentos construidos
    estacionario = modulo_raices(**params_coef)

    if estacionario:
        print("\nEl modelo es estacionario\n\n:)")
    else:
        print("\nEl modelo no es estacionario")
    


def invertible(modelo):
    "Esta función verifica si el modelo es invertible o no."
    ars = []
    elevados = []

    # Filtrar solo los parámetros que comienzan con 'ar'
    for i in modelo.params.index:
        if i.startswith('ma'):
            ars.append(modelo.params[i])
            grado = int(i.split('L')[-1])
            elevados.append(grado)

    # Construir los argumentos para modulo_raices_12.
    # Se asigna cada coeficiente a p{grado} y se dejan en cero los demás.
    params_coef = {}
    for coef, grado in zip(ars, elevados):
        params_coef[f"p{grado}"] = coef

    print("El polinomio de la parte media móvil es:", params_coef)
    

    print("")
    # Llamar a la función modulo_raices_12 con los argumentos construidos
    estacionario = modulo_raices(**params_coef)

    if estacionario:
        print("\nEl modelo es invertible\n\n:)")
    else:
        print("\nEl modelo no es invertible")
    


def parsimonia(modelo, alpha=0.10):
    "Esta función verifica si el modelo cumple el principio de parsimonia o no."
    if np.any(modelo.pvalues > alpha):
        print("Hay coeficientes no significativos, no se cumple el principio de parsimonia")
    else:
        print("El modelo cumple el principio de parsimonia\n\n:)")



# Función para calcular porcentaje dentro de k desviaciones estándar
def porcentaje_dentro_k(residuos, k, media, std):
    dentro = np.abs(residuos - media) <= k * std
    return np.mean(dentro) * 100

# Calculamos la media y desviación estándar
def normal_relajacion(residuos):
    # Cálculo para ±1σ, ±2σ, ±3σ
    for k, ref in zip([1, 2, 3], [68, 95, 99.7]):
        pct = porcentaje_dentro_k(residuos, k, np.mean(residuos), np.std(residuos))
        print(f"{pct:.2f}% de los residuos están dentro de ±{k}σ (esperado ≈ {ref}%)")

