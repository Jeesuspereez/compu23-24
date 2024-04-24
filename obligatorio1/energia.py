import matplotlib.pyplot as plt

def plot_from_file(filename, label):
    # Listas para almacenar los valores de x e y
    x_values = []
    y_values = []

    # Leer el archivo y procesar los datos
    with open(filename, 'r') as file:
        for line in file:
            # Eliminar espacios en blanco al inicio y al final de la línea
            line = line.strip()
            # Verificar si la línea está vacía
            if line:
                # Dividir la línea en dos partes usando la coma como separador
                parts = line.split(',')
                try:
                    # Convertir las partes en valores flotantes y agregarlos a las listas correspondientes
                    x_values.append(float(parts[0]))
                    y_values.append(float(parts[1]))
                except ValueError:
                    # Si hay un error al convertir a flotante, imprimir la línea que causó el error
                    print("Error en la línea:", line)

    # Graficar los datos
    plt.plot(x_values, y_values, label=label)
    
    # Devolver las listas de valores x e y
    return x_values, y_values

# Graficar desde cada archivo
x1, y1 = plot_from_file('potenciale.txt', 'Potencial')
x2, y2 = plot_from_file('cinetica.txt', 'Cinética')
x3, y3 = plot_from_file('energianplanet.txt', 'Energía')

# Etiquetas y título
plt.xlabel('Abscisa')
plt.ylabel('Energía')
plt.title('Gráfico de Energía vs Abscisa')

# Ajustar límites de los ejes x e y para que el origen coincida
min_x = min(min(x1), min(x2), min(x3))
max_x = max(max(x1), max(x2), max(x3))
min_y = min(min(y1), min(y2), min(y3))
max_y = max(max(y1), max(y2), max(y3))

# Ajustar límites para que incluyan el origen (0,0)
plt.xlim(left=min(0, min_x), right=max(0, max_x))
plt.ylim(bottom=min(0, min_y), top=max(0, max_y))

# Leyenda
plt.legend()

# Mostrar la gráfica
plt.grid(True)
plt.show()

