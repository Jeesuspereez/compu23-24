import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.signal import savgol_filter

# Leer datos de los archivos
#def leer_datos(nombre_archivo):
 #   with open(nombre_archivo, 'r') as file:
  #      data = [float(line.strip()) for line in file]
   # return data

   # Leer datos de los archivos
def leer_datos(nombre_archivo):
    with open(nombre_archivo, 'r') as file:
        data = [float(line.strip()) for line in file if line.strip()]
    return data

# Obtener datos de los archivos
tiempo = np.array(leer_datos('tiempo.txt'))
energia_cin = np.array(leer_datos('cineticalen.txt'))
energia_pot = np.array(leer_datos('potencialen.txt'))

# Ajustar los datos para que tengan la misma longitud
min_length = min(len(tiempo), len(energia_cin), len(energia_pot))
tiempo = tiempo[:min_length]
energia_cin = energia_cin[:min_length]
energia_pot = energia_pot[:min_length]

# Crear el gráfico
plt.figure(figsize=(18, 12))

# Interpolación de spline cúbico
spline_cin = CubicSpline(tiempo, energia_cin)
spline_pot = CubicSpline(tiempo, energia_pot)

# Evaluar las interpolaciones en un rango suave de tiempo
tiempo_smooth = np.linspace(tiempo.min(), tiempo.max(), 300)
energia_cin_smooth = spline_cin(tiempo_smooth)
energia_pot_smooth = spline_pot(tiempo_smooth)

# Calcular la suma de las dos curvas
energia_total_smooth = energia_cin_smooth + energia_pot_smooth

# Aplicar un suavizado adicional a la curva azul
energia_total_smooth = savgol_filter(energia_total_smooth, window_length=31, polyorder=3)

# Graficar las curvas resultantes
plt.plot(tiempo_smooth, energia_cin_smooth, label='Energía cinética', color='green')
plt.plot(tiempo_smooth, energia_pot_smooth, label='Energía potencial', color='red')

# Graficar la curva suavizada de la energía total
plt.plot(tiempo_smooth, energia_total_smooth, label='Energía total', color='blue')

# Etiquetas y leyenda
plt.xlabel('Tiempo', fontsize=25)
plt.ylabel('Energía', fontsize=25)
plt.title('Conservación de la energía', fontsize=25)
plt.legend(loc='best', fontsize=20)

# Ajustes estéticos
plt.grid(True)
plt.tick_params(axis="both", labelsize=20, labelrotation=0, labelcolor="black")

plt.show()
