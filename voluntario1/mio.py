import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt

PI = 3.141592

# Función para leer el archivo y extraer la tercera columna
def leer_tercera_columna(nombre_archivo):
    tercera_columna = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            columnas = linea.split()
            if len(columnas) >= 3:
                try:
                    valor = float(columnas[2])
                    tercera_columna.append(valor)
                except ValueError:
                    print(f"Error: No se pudo convertir a float en línea '{linea.strip()}'")
    return tercera_columna

# Función para leer la cuarta columna
def leer_cuarta_columna(nombre_archivo):
    cuarta_columna = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            columnas = linea.split()
            if len(columnas) >= 4:
                try:
                    valor = float(columnas[3])
                    cuarta_columna.append(valor)
                except ValueError:
                    print(f"Error: No se pudo convertir a float en línea '{linea.strip()}'")
    return cuarta_columna

# Función para leer la quinta columna
def leer_quinta_columna(nombre_archivo):
    quinta_columna = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            columnas = linea.split()
            if len(columnas) >= 5:
                try:
                    valor = float(columnas[4])
                    quinta_columna.append(valor)
                except ValueError:
                    print(f"Error: No se pudo convertir a float en línea '{linea.strip()}'")
    return quinta_columna

# Función de densidad de probabilidad P(v)
def P(v, m, Kb, T):
    return (m / (Kb * T)) * v * np.exp(-(m * v**2) / (2 * Kb * T))

# Función de densidad de probabilidad para vx
def Pejex(v, m, Kb, T):
    return sqrt(2 * m / (np.pi * Kb * T)) * np.exp(-(m * v**2) / (2 * Kb * T))

# Función para calcular y mostrar la densidad de probabilidad
def mostrar_densidad_probabilidad(datos1, datos2, datos3, m, Kb, T):
    plt.figure(figsize=(10, 6))

    # Histograma para la tercera columna
    plt.hist(datos1, bins=30, density=True, alpha=0.6, color='g', edgecolor='black', label='Densidad de Probabilidad módulo de velocidades')

    # Histograma para la cuarta columna
    plt.hist(datos2, bins=30, density=True, alpha=0.6, color='b', edgecolor='black', label='Densidad de Probabilidad velocidad eje x')

    # Histograma para la quinta columna
    plt.hist(datos3, bins=30, density=True, alpha=0.6, color='y', edgecolor='black', label='Densidad de Probabilidad velocidad eje y')

    # Calcular valores para P(v) y Pejex(v)
    x = np.linspace(min(min(datos1), min(datos2), min(datos3)), max(max(datos1), max(datos2), max(datos3)), 1000)
    y1 = P(x, m, Kb, T)
    y2 = Pejex(x, m, Kb, T)

    # Graficar P(v)
    plt.plot(x, y1, color='r', linewidth=2, label='Distribución del módulo de velocidades')

    # Graficar Pejex(v)
    plt.plot(x, y2, color='m', linewidth=2, label='Distribución de velocidad eje x y eje y')

    # Añadir etiquetas y título
    plt.title('Histograma del módulo de velocidades para T=0.99U (1U=120K)')
    plt.xlabel('Valores de la velocidad en UDS (1UD=1.58*10^2 m/s)')
    plt.ylabel('Densidad de probabilidad')
    plt.grid(True)
    plt.legend()
    plt.show()

# Parámetros para la función P(v)
m = 1.0    # masa (ejemplo)
Kb = 1.0  # constante de Boltzmann (ejemplo)
T = 0.9875631901  # temperatura (ejemplo)

# Nombre del archivo de texto
nombre_archivo = 'velocity.txt'

# Leer los datos de la tercera, cuarta y quinta columna
tercera_columna = leer_tercera_columna(nombre_archivo)
cuarta_columna = leer_cuarta_columna(nombre_archivo)
quinta_columna = leer_quinta_columna(nombre_archivo)

# Mostrar la densidad de probabilidad
mostrar_densidad_probabilidad(tercera_columna, cuarta_columna, quinta_columna, m, Kb, T)
