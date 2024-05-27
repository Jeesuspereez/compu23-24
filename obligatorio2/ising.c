//programa que simula metodo de ising para sistema magnetico
//hay que actualizar la matriz n^2 veces, entonces tenemos s1. Imprimimos s1. Volvemos hacer y tenemos s2, lo imprimimos. Hacer esto un millon de veces

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

//declaracion de funciones 
double aleatorio();
void matrizale(int **matriz, int n, int m); // Declaración de la función matrizale

int main(void)
{
    int **s, n, m, filas, columnas, pasomillon, i, j, nmenuno, mmenuno, nmasuno, mmasuno,N;
    double  energia, T, epsilon, p, expo;

    // Inicializar la semilla para generar números pseudoaleatorios
    srand(time(NULL));

     // Escribir la matriz en el archivo
    FILE *archivo;
    archivo = fopen("matriz.dat", "w");

    if (archivo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    //iniciamos contador pasomillon
    pasomillon=0;
    energia=0;

    // número de filas y columnas de la red OJO AL EMPEZAR EL PUNTERO EN CERO TENEMOS QUE HACER DIM-1
    filas=90;
    columnas=90;
    N=filas*filas;

    //probar a subir a 80x80 bajar fps y temp critica

    //Valor de la temperatura OJO T[0,5]
    T=4.;

    // Asignar memoria para las filas
    s = (int **)malloc((filas+1) * sizeof(int *));
    if (s == NULL) {
        printf("Error: no se pudo asignar memoria para las filas\n");
        return 1;
    }

    // Asignar memoria para las columnas de cada fila
    for (i = 0; i < columnas+1; i++) {
        s[i] = (int *)malloc((columnas+1) * sizeof(int));
        if (s[i] == NULL) {
            printf("Error: no se pudo asignar memoria para las columnas de la fila %d\n", i);
            return 1;
        }
    }

     // Inicializar la matriz (identidad)
  //  for (i = 0; i < filas; i++) {
   //    for (j = 0; j < columnas; j++) {
    //      s[i][j] = 1;
    //    }
  //  }

    //llamar a funcion matriz aleatoria
    matrizale(s, filas, columnas);

    //leer la matriz de fichero

     // Mostrar la matriz por pantalla
    printf("Matriz inicial:\n");
    for (i = 0; i < filas; i++) {
        for (j = 0; j < columnas; j++) {
            printf("%d ", s[i][j]);
        }
        printf("\n");
    }



//programa
for(int k=0;k<=10000; k++)
{
    //comenzamos un paso de montecarlo (N^2 iteraciones,en este caso N=filas*columnas)
    for(i=0; i<N ; i++)
    {
    //1.Generar posicion aleatoria en la red
    // Generar un número aleatorio para la fila (n) en el rango [0, filas - 1]
      n = rand() % filas;

    // Generar un número aleatorio para la columna (m) en el rango [0, columnas - 1]
      m = rand() % columnas;

    //2.Evaluar p
    //evaluar energia
   // energia= 2*s[n][m]*(s[n + 1][m] + s[n - 1][m] + s[n][m +1] + s[n][ m -1 ]);
   //formula 
    energia= 2*s[n][m]*(s[(n + 1)%filas][m] + s[(n - 1+filas)%filas][m] + s[n][(m +1)%columnas] + s[n][ (m -1+columnas)%columnas]);

        //evaluar exponencial
        expo=exp(-energia/T);

        //elegir p como el minimo entre 1 y la expo
        if(1<expo)
        {
            p=1;
        }
        else p=expo;

        //3.Crear numero aleatorio epsilon€[0,1], si eps<p cambia el signo del spin
  
       epsilon=aleatorio();
        if (epsilon<p)
        {
           s[n][m]= -s[n][m];
     
        }

         
    }

        //mostramos la matriz
      for (i = 0; i < filas; i++) {
        for (j = 0; j < columnas; j++) {
             if (j < columnas-1) {
            fprintf(archivo, "%d, ", s[i][j]);
             }
             else  fprintf(archivo, "%d ", s[i][j]);
        }
        fprintf(archivo, "\n");
    }

fprintf(archivo, "\n");

    }

    fclose(archivo);

    // Liberar la memoria
    for (i = 0; i < filas+1; i++) {
        free(s[i]);
    }
    free(s);

    return 0;
    }


//Implementacion de las funciones
  // Generar un número aleatorio en el rango [0, 1)
double aleatorio() 
{
    return (double)rand() / (double)RAND_MAX;
}

void matrizale(int **matriz, int n, int m) {
    int i, j, valor_aleatorio;
    
    for (i = 0; i < n; i++) { // Corregir el límite de iteración en función de las filas
        for (j = 0; j < m; j++) { // Corregir el límite de iteración en función de las columnas
            // Generar un número aleatorio entre 0 y 1 y luego convertirlo a +1 o -1
            valor_aleatorio = rand() % 2; // Genera 0 o 1
            matriz[i][j] = valor_aleatorio * 2 - 1; // Convierte 0 a -1 y 1 a +1
        }
    }
}