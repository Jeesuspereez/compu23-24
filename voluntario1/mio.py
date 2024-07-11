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

def leer_cuarta_columna(nombre_archivo):
    cuarta_columna = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            columnas = linea.split()
            if len(columnas) >= 4:  # Asegurarse de que haya al menos cuatro columnas en la línea
                cuarta_columna.append(float(columnas[3]))
    return cuarta_columna


# Función de densidad de probabilidad P(v)
def P(v, m, Kb, T):
    return (m / (Kb * T)) * v * np.exp(-(m * v**2) / (2 * Kb * T))

#Funcion de densidad de probabilidad para vx
def Pejex(v, m, Kb, T):
    return sqrt(2*m / (np.pi*Kb * T))  * np.exp(-(m * v**2) / (2 * Kb * T))


# Función para calcular y mostrar la densidad de probabilidad
def mostrar_densidad_probabilidad(datos, m, Kb, T):
    plt.figure(figsize=(10, 6))

    # Histograma
    plt.hist(datos, bins=30, density=True, alpha=0.6, color='g', edgecolor='black', label='Densidad de Probabilidad')

    # Calcular valores para P(v)
    x = np.linspace(min(datos), max(datos), 1000)
    y = P(x, m, Kb, T)
 #   y = Pejex(x, m, Kb, T)

    # Graficar P(v)
    plt.plot(x, y, color='r', linewidth=2, label='Distribución de Velocidades')

    # Añadir etiquetas y título
    plt.title('Histograma del módulo de velocidades para T=0.71U (1U=120K)')
  #  plt.title('Histograma del valor absoluto de la velocidad en el eje x para T=0.77U (1U=120K)')
    plt.xlabel('Valores de la velocidad en UDS (1UD=1.58*10^2 m/s)')
    plt.ylabel('Densidad de probabilidad')
    plt.grid(True)
    plt.legend()
    plt.show()

# Parámetros para la función P(v)
m = 1.0    # masa (ejemplo)
Kb = 1.0  # constante de Boltzmann (ejemplo)
T = 0.7075938680  # temperatura (ejemplo)

# Nombre del archivo de texto
nombre_archivo = 'velocity.txt'

# Leer los datos de la tercera columna
tercera_columna = leer_tercera_columna(nombre_archivo)
# cuarta_columna = leer_cuarta_columna(nombre_archivo)

# Mostrar la densidad de probabilidad
mostrar_densidad_probabilidad(tercera_columna, m, Kb, T)
# mostrar_densidad_probabilidad(cuarta_columna, m, Kb, T)
