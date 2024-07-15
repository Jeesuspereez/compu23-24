import matplotlib
import matplotlib.pyplot as plt;
import numpy as np;
import scipy
import scipy.optimize as opt;

plt.rcParams["font.family"] = "DejaVu Sans"
ax = plt.subplot()
# This is the function we are trying to fit to the data.
def func(x, a, b):
     return a*x+b

# Generate some data, you don't have to do this, as you already have your data
archivodatos = np.loadtxt('TyP.txt')
xdata=archivodatos[:,0]
ydata=archivodatos[:,1]

x = np.linspace(min(xdata), max(xdata), 100)

# Plot the actual data
ax.set_xlabel('Temperatura', fontname='DejaVu Sans', fontsize='12')
ax.set_ylabel('Presión', fontname='DejaVu Sans', fontsize='12')
ax.plot(xdata, ydata, ".", label="Data", color='#3A405A', markersize=10)

# The actual curve fitting happens here
optimizedParameters, pcov = opt.curve_fit(func, xdata, ydata)

print(*optimizedParameters)
print(np.sqrt(np.diag(pcov)))

# Use the optimized parameters to plot the best fit
ax.plot(x, func(x, *optimizedParameters), label="fit", color='#B4045F')

# Show the graph
plt.savefig('Ajuste.pdf',dpi=300, bbox_inches = "tight")