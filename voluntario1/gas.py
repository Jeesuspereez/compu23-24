import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

# Función para leer los datos del archivo y separarlos en instantes de tiempo
def leer_datos(archivo):
    datos = []
    instante = []
    with open(archivo, 'r') as f:
        for linea in f:
            if linea.strip():  # Si la línea no está en blanco
                x, y = map(float, linea.split(','))
                instante.append((x, y))
            else:  # Si encontramos un salto de línea, añadimos los datos al resultado y reiniciamos para el siguiente instante
                datos.append(instante)
                instante = []
    return datos

# Leer datos del archivo
datos = leer_datos('lennard.txt')

# Configuración de la animación
fig, ax = plt.subplots()
sc = ax.scatter([], [])
ax.set_xlim(0, 10)  # Ajusta los límites x e y según tus datos
ax.set_ylim(0, 10)
ax.set_xlabel('Posición x')
ax.set_ylabel('Posición y')
ax.set_title('Movimiento de las partículas')

# Función de inicialización de la animación
def init():
    sc.set_offsets(np.array([]).reshape(0, 2))  # Inicializar con un array vacío
    return sc,

# Función de actualización de la animación
def update(frame):
    if frame < len(datos):
        x = [particula[0] for particula in datos[frame]]
        y = [particula[1] for particula in datos[frame]]
        sc.set_offsets(np.column_stack((x, y)))
    return sc,

# Crear la animación
ani = FuncAnimation(fig, update, frames=len(datos), init_func=init, blit=True)

# Opción para guardar el video o mostrarlo en pantalla
opcion = input("¿Deseas guardar el video (s/n)? ").lower()
if opcion == 's':
    nombre_video = input("Ingresa el nombre del archivo de video: ")
    ani.save(nombre_video, writer='ffmpeg')
    print(f"Video guardado como {nombre_video}")
else:
    plt.show()