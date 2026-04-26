import numpy as np
from scipy.stats import ttest_1samp

def media_cero(modelo):
    residuos=modelo.resid
    pvalor = ttest_1samp(residuos, 0)[1] # Guardamos el p-valor
    if pvalor < 0.05:
        print(f"Los residuos no son cero, pvalor={pvalor}, media={residuos.mean()}")
    else:
        print(f"Los residuos son cero, pvalor={pvalor}, media={residuos.mean()}\n\n:)")

def modulo_raices_48(
    p1=0, p2=0, p3=0, p4=0, p5=0, p6=0, p7=0, p8=0, p9=0, p10=0, 
    p11=0, p12=0, p13=0, p14=0, p15=0, p16=0, p17=0, p18=0, p19=0, p20=0,
    p21=0, p22=0, p23=0, p24=0, p25=0, p26=0, p27=0, p28=0, p29=0, p30=0,
    p31=0, p32=0, p33=0, p34=0, p35=0, p36=0, p37=0, p38=0, p39=0, p40=0,
    p41=0, p42=0, p43=0, p44=0, p45=0, p46=0, p47=0, p48=0
):
    coef = np.zeros(49)

    # Asignar los coeficientes en orden inverso
    coef[0] = p48  # x^48
    coef[1] = p47  # x^47
    coef[2] = p46  # x^46
    coef[3] = p45  # x^45
    coef[4] = p44  # x^44
    coef[5] = p43  # x^43
    coef[6] = p42  # x^42
    coef[7] = p41  # x^41
    coef[8] = p40  # x^40
    coef[9] = p39  # x^39
    coef[10] = p38 # x^38
    coef[11] = p37 # x^37
    coef[12] = p36 # x^36
    coef[13] = p35 # x^35
    coef[14] = p34 # x^34
    coef[15] = p33 # x^33
    coef[16] = p32 # x^32
    coef[17] = p31 # x^31
    coef[18] = p30 # x^30
    coef[19] = p29 # x^29
    coef[20] = p28 # x^28
    coef[21] = p27 # x^27
    coef[22] = p26 # x^26
    coef[23] = p25 # x^25
    coef[24] = p24 # x^24
    coef[25] = p23 # x^23
    coef[26] = p22 # x^22
    coef[27] = p21 # x^21
    coef[28] = p20 # x^20
    coef[29] = p19 # x^19
    coef[30] = p18 # x^18
    coef[31] = p17 # x^17
    coef[32] = p16 # x^16
    coef[33] = p15 # x^15
    coef[34] = p14 # x^14
    coef[35] = p13 # x^13
    coef[36] = p12 # x^12
    coef[37] = p11 # x^11
    coef[38] = p10 # x^10
    coef[39] = p9  # x^9
    coef[40] = p8  # x^8
    coef[41] = p7  # x^7
    coef[42] = p6  # x^6
    coef[43] = p5  # x^5
    coef[44] = p4  # x^4
    coef[45] = p3  # x^3
    coef[46] = p2  # x^2
    coef[47] = p1  # x^1
    coef[48] = 1   # x^0

    raices = np.roots(coef)  # Raíces del polinomio
    modulo = np.abs(raices)  # Módulo de las raíces
    estacionario = np.all(modulo > 1)  # Verificar si todas las raíces están fuera del círculo unitario

    print("Raíces del polinomio característico:", raices)
    print("\nMódulo de las raíces:", modulo)
    print("\n¿Las raíces están fuera del círculo unitario? ", estacionario)

    return estacionario



def estacionario(modelo):

    ars = []
    elevados = []

    # Filtrar solo los parámetros que comienzan con 'ar'
    for name, coef in zip(modelo.param_names, modelo.params):
        if name.startswith('ar'):
            ars.append(-coef)
            grado = int(name.split('L')[-1])
            elevados.append(grado)

    # Construir los argumentos para modulo_raices_12.
    # Se asigna cada coeficiente a p{grado} y se dejan en cero los demás.
    params_coef = {}
    for coef, grado in zip(ars, elevados):
        params_coef[f"p{grado}"] = coef

    print("El polinomio de la parte autoregresiva es:", params_coef)
    print("")
    
    # Llamar a la función modulo_raices_48 con los argumentos construidos
    estacionario = modulo_raices_48(**params_coef)

    if estacionario:
        print("\nEl modelo es estacionario\n\n:)")
    else:
        print("\nEl modelo no es estacionario")





def invertible(modelo):

    ars = []
    elevados = []

    # Filtrar solo los parámetros que comienzan con 'ar'
    for name, coef in zip(modelo.param_names, modelo.params):
        if name.startswith('ma'):
            ars.append(coef)
            grado = int(name.split('L')[-1])
            elevados.append(grado)

    # Construir los argumentos para modulo_raices_12.
    # Se asigna cada coeficiente a p{grado} y se dejan en cero los demás.
    params_coef = {}
    for coef, grado in zip(ars, elevados):
        params_coef[f"p{grado}"] = coef

    print("El polinomio de la parte media móvil es:", params_coef)
    print("")
    
    # Llamar a la función modulo_raices_48 con los argumentos construidos
    estacionario = modulo_raices_48(**params_coef)

    if estacionario:
        print("\nEl modelo es invertible\n\n:)")
    else:
        print("\nEl modelo no es invertible")
    



def parsimonia(modelo, alpha=0.10):
    if np.any(modelo.pvalues > alpha):
        print("Hay coeficientes no significativos, no se cumple el principio de parsimonia")
    else:
        print("El modelo cumple el principio de parsimonia\n\n:)")