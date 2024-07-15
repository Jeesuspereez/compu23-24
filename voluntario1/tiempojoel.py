import subprocess
import time

# Guarda la hora de inicio
inicio = time.time()

# Ejecuta el programa usando subprocess
subprocess.run(["./nplanet.exe"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

# Guarda la hora de finalización
fin = time.time()

# Calcula la diferencia de tiempo
tiempo_transcurrido = fin - inicio

# Imprime el tiempo de ejecución
print("Tiempo de ejecución:", tiempo_transcurrido, "segundos")
