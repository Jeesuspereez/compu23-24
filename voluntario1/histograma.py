import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell

# Leer datos de velocidades
velocities = np.loadtxt('velocity.txt')
v_magnitudes = np.sqrt(velocities[:, 0]**2 + velocities[:, 1]**2)

# Calcular la temperatura del archivo de energía (opcional)
#temperature = 0.5 * np.mean(v_magnitudes**2)
temperature=0.0170648522

# Parámetros de la distribución de Maxwell-Boltzmann
kB = 1.0
m = 1.0
scale = np.sqrt(kB * temperature / m)

# Histograma de las velocidades
#plt.hist(v_magnitudes, bins=30, density=True, alpha=0.6, color='g', label='Simulación')

# Distribución de Maxwell-Boltzmann teórica
x = np.linspace(0, np.max(v_magnitudes), 100)
pdf = maxwell.pdf(x, scale=scale)
plt.plot(x, pdf, 'r-', label='Maxwell-Boltzmann teórica')

plt.xlabel('Velocidad')
plt.ylabel('Densidad de probabilidad')
plt.title('Comparación del histograma de velocidades con la distribución de Maxwell-Boltzmann')
plt.legend()
plt.show()
