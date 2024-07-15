import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# Datos
V = np.array([0.009869604401, 0.0197392088, 0.0296088132, 0.04934802201, 0.09869604401, 0.4934802201, 0.9869604401])
K = np.array([0.939, 0.934, 0.9216666667, 0.8276666667, 0.061, 0, 0])
Error_K = np.array([0.005507570547, 0.001, 0.00202758751, 0.005174724899, 0.003511884584, 0, 0])

# Interpolación lineal
interpolacion_lineal = interp1d(V, K, kind='linear')

# Creando nuevos puntos para la interpolación
V_nuevos = np.linspace(V.min(), V.max(), 100)
K_interpolados = interpolacion_lineal(V_nuevos)

# Graficar los datos y la interpolación
plt.errorbar(V, K, yerr=Error_K, fmt='o', label='Datos con error', capsize=5)
plt.plot(V_nuevos, K_interpolados, '-', label='Interpolación lineal')
plt.xlabel('V(x)')
plt.ylabel('K')
plt.title('K frente a V con N=1000 y nciclos=50')
plt.legend()
plt.grid(True)
plt.show()
