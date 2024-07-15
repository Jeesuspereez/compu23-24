import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# Datos de la tabla
N = np.array([50000, 100000, 150000, 200000, 250000, 300000])
pc = np.array([4.968, 9.827, 14.5, 19.9, 24.937, 30.168])
O2 = np.array([4.9, 9.742, 14.42, 19.5, 24.368, 29.012])
O3 = np.array([4.89, 9.73, 14.418, 19.1, 24.1, 28.0])

# Crear spline para pc
spl_pc = make_interp_spline(N, pc)
N_smooth_pc = np.linspace(N.min(), N.max(), 500)
pc_smooth = spl_pc(N_smooth_pc)

# Crear spline para O2
spl_O2 = make_interp_spline(N, O2)
N_smooth_O2 = np.linspace(N.min(), N.max(), 500)
O2_smooth = spl_O2(N_smooth_O2)

# Crear spline para O3
spl_O3 = make_interp_spline(N, O3)
N_smooth_O3 = np.linspace(N.min(), N.max(), 500)
O3_smooth = spl_O3(N_smooth_O3)

# Graficar las funciones
plt.figure(figsize=(10, 6))
plt.plot(N_smooth_pc, pc_smooth, label='Pc', color='blue')
plt.plot(N_smooth_O2, O2_smooth, label='O2', color='green')
plt.plot(N_smooth_O3, O3_smooth, label='O3', color='red')
plt.scatter(N, pc, color='blue', marker='o')
plt.scatter(N, O2, color='green', marker='x')
plt.scatter(N, O3, color='red', marker='s')
plt.xlabel('P iteraciones')
plt.ylabel('Tiempo')
plt.title('Comparación de pc, O2 y O3 para p iteraciones')
plt.legend()
plt.grid(True)
plt.show()
