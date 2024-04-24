import matplotlib.pyplot as plt

def plot_from_file(filename, label):
    x_values = []
    y_values = []
    with open(filename, 'r') as file:
        for line in file:
            # Verificar si la línea está vacía o contiene solo espacios en blanco
            if line.strip():
                parts = line.split(',')
                x_values.append(float(parts[0]))
                y_values.append(float(parts[1]))
    plt.plot(x_values, y_values, label=label)

# Graficar desde cada archivo
plot_from_file('energianplanet.txt', 'Energía')
plot_from_file('cinetica.txt', 'Cinética')
plot_from_file('potenciale.txt', 'Potencial')

# Etiquetas y título
plt.xlabel('Abscisa')
plt.ylabel('Energía')
plt.title('Gráfico de Energía vs Abscisa')

# Ajustar límites de los ejes x e y
all_x_values = []
all_y_values = []
for filename in ['energianplanet.txt', 'cinetica.txt', 'potenciale.txt']:
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():
                parts = line.split(',')
                all_x_values.append(float(parts[0]))
                all_y_values.append(float(parts[1]))

# Calcular el rango de los valores de los ejes x e y
x_range = max(all_x_values) - min(all_x_values)
y_range = max(all_y_values) - min(all_y_values)

# Calcular el rango máximo entre los ejes x e y
max_range = max(x_range, y_range)

# Establecer los límites de los ejes x e y centrados en cero
plt.xlim(0, max_range*1000000)
plt.ylim(-max_range*1000000/2, 1000000*max_range/2)

# Ajustar la posición del eje x en la mitad de la figura
plt.gca().spines['bottom'].set_position(('data', 0))

# Mostrar la leyenda
plt.legend()

# Mostrar la gráfica
plt.grid(True)
plt.show()


