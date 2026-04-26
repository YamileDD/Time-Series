import numpy as np

# Definir la matriz de autocorrelación R
R = np.array([
    [1, 0.2, 0.5],
    [0.2, 1, 0.2],
    [0.5, 0.2, 1]
])

# Definir el vector de autocorrelaciones rho
rho = np.array([0.2, 0.5, 0.7])

# Resolver el sistema de ecuaciones para los parámetros phi_1, phi_2 y phi_3
phi = np.linalg.solve(R, rho)

print(phi)
