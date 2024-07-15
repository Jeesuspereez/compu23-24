import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Datos de la tabla
N = np.array([10, 20, 30, 40, 50, 60, 70, 80, 90, 100])
Pc = np.array([7.248, 13.269, 20.762, 28.908, 39.439, 51.946, 66.06, 83.616, 102.033, 117.155])
Joel = np.array([9.64, 18.96, 30.92, 46.48, 65.31, 86.45, 111.11, 138.47, 170.14, 205.71])

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
plt.plot(N_smooth_Pc, Pc_smooth, label='Pc AMD Ryzen 7 6800HS', color='blue')
plt.plot(N_smooth_Joel, Joel_smooth, label='Joel', color='green')
plt.scatter(N, Pc, color='blue', marker='o')
plt.scatter(N, Joel, color='green', marker='x')
plt.xlabel('N')
plt.ylabel('Tiempo (s)')
plt.title('Gráfica de las funciones Pc y Joel para sistema solar de N Jupiter')
plt.legend()
plt.grid(True)
plt.show()
