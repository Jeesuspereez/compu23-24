import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Datos de la tabla
N = np.array([5, 25, 35, 40, 45, 50, 60])
Joel = np.array([1.261, 10.148, 18.793, 24.343, 31.301, 38.304, 57.931])
Pc = np.array([1.72, 12.391, 22.936, 30.485, 40.368, 51.091, 82.121])

# Crear spline para Pc
spl_Pc = make_interp_spline(N, Pc)
N_smooth_Pc = np.linspace(N.min(), N.max(), 500)
Pc_smooth = spl_Pc(N_smooth_Pc)

# Crear spline para Joel
spl_Joel = make_interp_spline(N, Joel)
N_smooth_Joel = np.linspace(N.min(), N.max(), 500)
Joel_smooth = spl_Joel(N_smooth_Joel)

# Graficar las funciones
plt.figure(figsize=(10, 6))
plt.plot(N_smooth_Pc, Pc_smooth, label='Pc', color='blue')
plt.plot(N_smooth_Joel, Joel_smooth, label='Joel', color='green')
plt.scatter(N, Pc, color='blue', marker='o')
plt.scatter(N, Joel, color='green', marker='x')
plt.xlabel('N')
plt.ylabel('Tiempo')
plt.title('Comparación de Joel y Pc para N particulas en el gas')
plt.legend()
plt.grid(True)
plt.show()
