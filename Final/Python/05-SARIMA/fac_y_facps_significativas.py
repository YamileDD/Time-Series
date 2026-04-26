import numpy as np

def significancia(N, rhos):
    return 2 * np.sqrt(1 / (N) * (1 + 2 * sum(rho**2 for rho in rhos)))


def FAC(N: int, rhos: list, print_rhos: bool = False) -> list:
    """
    Nos da la lista de los rhos significativos
    """
    significativos = []
    rho_list = [0] * len(rhos)
    
    updated_rhos = list(rho_list)  # Copy to retain original values when updating
    for i, r in enumerate(rhos):  # for index, rho in enumerate(r_list):
        threshold_value = significancia(N, updated_rhos)    # se actualiza el valor con los rhos significativos
        if abs(r) > threshold_value:
            if print_rhos: 
                print(f"rho {i+1} es significativo")
                print(f"{abs(r)} > {threshold_value}\n")

            updated_rhos[i] = r
            significativos.append((r, i+1))
        else:
            if print_rhos: 
                print(f"rho {i+1} no es significativo")
                print(f"{abs(r)} < {threshold_value}\n")

    print("Valores de autocorrelacion significativos:")
    for val, i in significativos:
        print(f"r{i}: {val}")

    return updated_rhos, significativos

def intervalo(N):
    "Esta es de las autocorrelaciones parciales"
    return 2*np.sqrt(1/ N)

def FACP(N: int, rho_list: list, print_rhos: bool = False) -> list:
    """
    Esta función nos da la lista de los rhos significativos
    """
    
    fueradelintervalo = intervalo(N)
    significativos = []

    for i, r in enumerate(rho_list):
        
        if abs(r) > fueradelintervalo:
            if print_rhos: 
                print(f"rho {i+1} es significativo")
                print(f"{abs(r)} > {fueradelintervalo}\n")
            significativos.append((r, i+1))
        else:
            if print_rhos: 
                print(f"rho {i+1} no es significativo")
                print(f"{abs(r)} < {fueradelintervalo}\n")
                pass

    print(f"Valores de autocorrelacion parcial significativos:")
    for r, i in significativos:
        print(f"rho {i}: {r}")

    return significativos