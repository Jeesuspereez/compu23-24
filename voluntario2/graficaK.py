import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Datos
N = np.array([500, 1000, 2000, 3000, 4000, 5000, 6000])
K = np.array([0.9223333333, 0.833, 0.5463333333, 0.3743333333, 0.279, 0.2403333333, 0.2406666667])
Error_K = np.array([0.003382963855, 0.004358898944, 0.006887992773, 0.01074450764, 0.006244997998, 0.003711842909, 0.0197681])

# Interpolación lineal
interpolacion_lineal = interp1d(N, K, kind='linear')

# Creando nuevos puntos para la interpolación
N_nuevos = np.linspace(500, 6000, 100)
K_interpolados = interpolacion_lineal(N_nuevos)

# Graficar los datos y la interpolación
plt.errorbar(N, K, yerr=Error_K, fmt='o', label='Datos con error', capsize=5)
plt.plot(N_nuevos, K_interpolados, '-', label='Interpolación lineal')
plt.xlabel('N')
plt.ylabel('K')
plt.title('Interpolación de K frente a N')
plt.legend()
plt.grid(True)
plt.show()
