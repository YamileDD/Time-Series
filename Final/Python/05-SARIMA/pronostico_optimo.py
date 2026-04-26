from scipy.stats import norm
import numpy as np

def po(H, modelo, ypre, res, t0, alpha=0.05):
    """
    Genera pronósticos recursivos con intervalos de confianza normales.
    
    Parámetros:
    - H: horizonte de pronóstico
    - modelo: modelo SARIMA ajustado (statsmodels)
    - ypre: lista de valores observados hasta t0
    - res: residuos históricos hasta t0
    - t0: índice del último valor observado
    - alpha: nivel de significancia (0.05 → IC del 95%)

    Retorna:
    - Listas: pronósticos, límites inferiores, límites superiores
    """
    # Estimación de varianza del error
    sigma = np.std(res, ddof=1)
    z = norm.ppf(1 - alpha / 2)

    # Inicializar listas
    y_forecast = list(ypre)
    res_forecast = list(res)
    forecasts, lower_bounds, upper_bounds = [], [], []

    # Parámetros del modelo
    phi1 = modelo.params[0]
    phi2 = modelo.params[1]
    theta1 = modelo.params[2]
    theta2 = modelo.params[3]
    theta3 = modelo.params[4]
    phi_s = modelo.params[5]
    theta_s = modelo.params[6]

    for h in range(1, H + 1):
        t = t0 + h

        # Función de respaldo para índices negativos
        def safe_get(lst, idx):
            return lst[idx] if idx >= 0 else 0.0

        # Pronóstico determinista
        pred = (
            phi1 * safe_get(y_forecast, t - 1) +
            phi2 * safe_get(y_forecast, t - 2) +
            phi_s * safe_get(y_forecast, t - 12) +
            theta1 * safe_get(res_forecast, t - 1) +
            theta2 * safe_get(res_forecast, t - 2) +
            theta3 * safe_get(res_forecast, t - 3) +
            theta_s * safe_get(res_forecast, t - 12)
        )

        # Varianza aproximada del pronóstico (simple: sigma^2 * h)
        var_h = sigma**2 * h
        margin = z * np.sqrt(var_h)

        # Guardar resultados
        forecasts.append(pred)
        lower_bounds.append(pred - margin)
        upper_bounds.append(pred + margin)

        y_forecast.append(pred)
        res_forecast.append(0.0)  # Suponemos errores futuros en media cero

    return forecasts, lower_bounds, upper_bounds
