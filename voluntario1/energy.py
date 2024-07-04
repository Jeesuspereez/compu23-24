import matplotlib.pyplot as plt

def leer_datos(archivo):
    with open(archivo, 'r') as f:
        datos = f.read().splitlines()
    datos = [float(i) for i in datos]
    return datos

# Leer datos de los archivos
tiempo = leer_datos('tiempo.txt')
cinetica = leer_datos('cineticalen.txt')
potencial = leer_datos('potencialen.txt')
energia = leer_datos('energianlen.txt')

# Crear la gráfica
plt.figure(figsize=(10, 6))

plt.plot(tiempo, cinetica, label='Energía Cinética', color='r')
plt.plot(tiempo, potencial, label='Energía Potencial', color='g')
plt.plot(tiempo, energia, label='Energía Total', color='b')

plt.xlabel('Tiempo')
plt.ylabel('Energía')
plt.title('Energías en función del tiempo')
plt.legend()

# Mostrar la gráfica
plt.show()
