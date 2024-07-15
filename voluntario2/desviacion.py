import numpy as np

def calcular_media_y_desviacion_estandar(file_path):
    # Leer los datos del archivo de texto
    data = np.loadtxt(file_path)
    
    # Calcular la media
    media = np.mean(data)
    
    # Calcular la desviación estándar
    desviacion_estandar = np.std(data)
    
    # Imprimir la media y la desviación estándar
    print(f'Media: {media}')
    print(f'Desviación Estándar: {desviacion_estandar}')

# Ruta al archivo de texto
file_path = 'cinetica.txt'

# Llamar a la función
calcular_media_y_desviacion_estandar(file_path)
