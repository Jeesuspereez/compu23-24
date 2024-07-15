import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import pearsonr

# Datos de la tabla
barreras = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])
k = np.array([0.9533333333, 0.893, 0.8346666667, 0.6856666667, 0.671, 0.492, 0.4546666667, 0.4, 0.2766666667, 0.1603333333, 0.04666666667])
error_k = np.array([0.004096068576, 0.008736894948, 0.002728450924, 0.001201850425, 0.004041451884, 0.003214550254, 0.002185812841, 0.004582575695, 0.00348010217, 0.002185812841, 0.0202758751])

# Función exponencial para el ajuste
def func_exp(x, a, b):
    return a * np.exp(b * x)

# Ajuste por mínimos cuadrados
popt, pcov = curve_fit(func_exp, barreras, k)
a, b = popt

# Valores ajustados
k_fit = func_exp(barreras, a, b)

# Coeficiente de correlación de Pearson
pearson_corr, _ = pearsonr(k, k_fit)

# Gráfica de los datos y el ajuste
plt.errorbar(barreras, k, yerr=error_k, fmt='o', label='Datos experimentales', capsize=5)
plt.plot(barreras, k_fit, label=f'Ajuste exponencial: $k = {a:.3f}e^{{{b:.3f}x}}$', color='red')
plt.xlabel('Barreras')
plt.ylabel('K')
plt.title('Ajuste Exponencial')
plt.legend()
plt.grid(True)
plt.show()

# Imprimir los resultados
print(f'Ajuste exponencial: k = {a:.3f} * exp({b:.3f} * x)')
print(f'Coeficiente de correlación de Pearson: {pearson_corr:.3f}')
