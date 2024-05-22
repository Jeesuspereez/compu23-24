import matplotlib.pyplot as plt
import numpy as np

# Leer los datos del archivo norma.dat
with open('norma.dat', 'r') as file:
    data = file.readlines()

# Convertir los datos a float
values = [float(line.strip()) for line in data]

# Crear el eje N
N = np.arange(1, len(values) + 1)

# Graficar los datos
plt.figure(figsize=(10, 6))
plt.plot(N, values, marker='o', linestyle='-', color='b', label='Valores de norma.dat')
plt.xlabel('Cada iteración')
plt.ylabel('Norma')
plt.title('Conservación de la norma en función del tiempo')
plt.legend()
plt.grid(True)
plt.show()
