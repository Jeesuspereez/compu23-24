import matplotlib.pyplot as plt

def leer_archivo(nombre_archivo):
    puntos_x = []
    puntos_y = []
    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            x, y = linea.strip().split(',')
            puntos_x.append(float(x))
            puntos_y.append(float(y))
    return puntos_x, puntos_y

def graficar_puntos(puntos_x, puntos_y):
    plt.scatter(puntos_x, puntos_y, color='blue', label='Puntos')
    plt.xlabel('Abscisa (x)')
    plt.ylabel('Coordenada (y)')
    plt.title('Gráfico de puntos')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    nombre_archivo = input("Ingrese el nombre del archivo .txt: ")
    puntos_x, puntos_y = leer_archivo(nombre_archivo)
    graficar_puntos(puntos_x, puntos_y)

