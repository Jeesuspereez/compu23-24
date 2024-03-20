//programa simulacion de planetas

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

#define G 6,67E-11



// Función para sumar h a cada elemento de una matriz de n filas y 2 columnas
void posicion(double matriz[][2], int n, double h) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2; j++) {
            // Sumar h a cada elemento de la matriz
            matriz[i][j] += h;
        }
    }
}

int main() {
    // Ejemplo de matriz de n filas y 2 columnas
       int n = 5;
    double matriz[][2] = {{1.0, 2.0},
                          {3.0, 4.0},
                          {5.0, 6.0},
                          {7.0, 8.0},
                          {9.0, 10.0}};
    double h = 2.0;

    // Mostrar la matriz antes de la suma
    printf("Matriz antes de sumar h:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%.2f ", matriz[i][j]);
        }
        printf("\n");
    }

    // Llamar a la función para sumar h a cada elemento de la matriz
    posicion(matriz, n, h);

    // Mostrar la matriz después de la suma
    printf("Matriz después de sumar h:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%.2f ", matriz[i][j]);
        }
        printf("\n");
    }

    return 0;
}
