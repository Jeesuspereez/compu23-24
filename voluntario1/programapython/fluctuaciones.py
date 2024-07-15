import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

def leer_datos(nombre_archivo):
    with open(nombre_archivo, 'r') as file:
        data = [float(line.strip()) for line in file if line.strip()]
    return data

# Función para seleccionar archivo
def seleccionar_archivo():
    print("Seleccione el archivo a visualizar:")
    print("1. Fluctuaciones")
    print("2. Temperatura")
    print("3. Correlación")
    opcion = input("Ingrese el número de su opción: ")
    if opcion == '1':
        return 'fluctuaciones.txt'
    elif opcion == '2':
        return 'temperatura.txt'
    elif opcion == '3':
        return 'correlacion.txt'
    else:
        print("Opción no válida. Mostrando fluctuaciones por defecto.")
        return 'fluctuaciones.txt'

# Seleccionar archivo
archivo_seleccionado = seleccionar_archivo()

# Obtener datos de los archivos
tiempo = np.array(leer_datos('tiempo.txt'))
energia = np.array(leer_datos(archivo_seleccionado))

# Ajustar los datos para que tengan la misma longitud
min_length = min(len(tiempo), len(energia))
tiempo = tiempo[:min_length]
energia = energia[:min_length]

# Crear el gráfico
plt.figure(figsize=(18, 12))

# Interpolación de spline cúbico
spline = CubicSpline(tiempo, energia)

# Evaluar las interpolaciones en un rango suave de tiempo
tiempo_smooth = np.linspace(tiempo.min(), tiempo.max(), 300)
energia_smooth = spline(tiempo_smooth)

# Graficar la curva interpolada
plt.plot(tiempo_smooth, energia_smooth, label='Temperatura', color='red')

# Etiquetas y leyenda
plt.xlabel('Tiempo', fontsize=25)
plt.ylabel('Temperatura', fontsize=25)
plt.title('Temperatura al calentar un factor de 1.1 en t=180 ', fontsize=25)
plt.legend(loc='best', fontsize=20)

# Ajustes estéticos
plt.grid(True)
plt.tick_params(axis="both", labelsize=20, labelrotation=0, labelcolor="black")

plt.show()
