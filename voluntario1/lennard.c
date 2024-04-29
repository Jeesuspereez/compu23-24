/* programa que simula la dinámica molecular de un gas con un potencial de lennard jones
-V(r)=4eps[(sigma/r)^12 - (sigma/r)^6]
-F(r)= -parcialV/r = 24eps*sigma^6(2sigma^6-r^6)/r^13 */  

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define Argonmass (6.6335209 * pow(10, -26)) // Masa en kg
#define Armstrong (1 * pow(10, -10)) // Armstrong en metros

// Declaración de funciones
void rescm(double *masa, int n);
void rescr(double *pos, int sig, int n);
void corriger(double *pos, int sig, int n);
void aceleracion(double *aceleracion, double *posx, double *posy, int L, int n);
void calculopos(double *pos, double *vel, double *acel, double h, int n);
void calculov(double *velocidad, double *aux, double *acel, double h, int n);
void calculow(double *aux, double *vel, double *acel, double h, int n);
double cinetica(double *masa, double *velx, double *vely, int n);
double potencial(double *masa, double *posx, double *posy, double *acelx, double *acely, int n);
double mangular(double *masa, double *posx, double *posy, double *velx, double *vely, int n);
void generar_posiciones(double *pos_x, double *pos_y, int dimension, int tamano_red);
void generar_velocidades(double *velx, double *vely, int dimension);
void imprimirCoordenadas(double x[], double y[], int n);

int main(void)
{
    // Abriendo ficheros
    FILE *archivo;
    char simplanet[] = "lennard.txt";
    archivo = fopen(simplanet, "w");

    FILE *archivo_cinetic;
    char cinetic[] = "cineticalen.txt";
    archivo_cinetic = fopen(cinetic, "w");

    FILE *archivo_potenciale;
    char potenciale[] = "potencialen.txt";
    archivo_potenciale = fopen(potenciale, "w");

    FILE *archivo_;
    char conservacion[] = "energianlen.txt";
    archivo_ = fopen(conservacion, "w");

    FILE *archivo__;
    char momentoangular[] = "angularlen.txt";
    archivo__ = fopen(momentoangular, "w");

    // Definición de variables
    int i, filas, j, k, N, L, sigma, epsilon;
    double *m, *r_x, *r_y, *v_x, *w_x, *w_y, *v_y, *a_x, *a_y, t, h;

    t = 0;
    h=0.001;
    sigma = 1;
    epsilon = 1;
    L = 10;
    N = 10;
    filas = N;

    // Asignando memoria para los vectores
    m = (double *)malloc(filas * sizeof(double));
    r_x = (double *)malloc(filas * sizeof(double));
    r_y = (double *)malloc(filas * sizeof(double));
    v_x = (double *)malloc(filas * sizeof(double));
    v_y = (double *)malloc(filas * sizeof(double));
    w_x = (double *)malloc(filas * sizeof(double));
    w_y = (double *)malloc(filas * sizeof(double));
    a_x = (double *)malloc(filas * sizeof(double));
    a_y = (double *)malloc(filas * sizeof(double));

    // Inicializando aceleraciones a cero
    for (i = 0; i < filas; i++)
    {
        r_x[i] = 0;
        r_y[i] = 0;
        v_x[i] = 0;
        v_y[i] = 0;
        a_x[i] = 0;
        a_y[i] = 0;
        
    }

    //asignamos posicion random
     generar_posiciones(r_x, r_y, N, L);

    // Rescalando vectores de posición
    rescr(r_x, sigma, filas);
    rescr(r_y, sigma, filas);

    // Inicializando masas
    for (i = 0; i < N; i++)
    {
        m[i] = Argonmass;
    }

    // Rescalamos masas
    rescm(m, filas);

    generar_velocidades(v_x,v_y,filas);

    // Rescalando velocidades
   /* rescr(v_y, sigma, filas);
    corriget(v_y, filas);*/

    imprimirCoordenadas(r_x, r_y, filas);

    // Calculando aceleraciones para la primera iteración
    aceleracion(a_x, r_x, r_y, L, filas);
    aceleracion(a_y, r_y, r_x, L, filas); 

    //mostramos vectores posiciones por pantalla a ver si esta bien


    /* ALGORITMO DE VERLET */
    for (t = 0; t < 100; t++)
    {

        // Calculando nuevas posiciones
        calculopos(r_x, v_x, a_x, t, filas);
        calculopos(r_y, v_y, a_y, t, filas);

        // Calculando w
        calculow(w_x, v_x, a_x, t, filas);
        calculow(w_y, v_y, a_y, t, filas);

        // Calculando nuevas aceleraciones
        aceleracion(a_x, r_x, r_y, L, filas);
        aceleracion(a_y, r_y, r_x, L, filas);

        // Calculando nuevas velocidades
        calculov(v_x , w_x, a_x, t, filas);
        calculov(v_y , w_y, a_y, t, filas);

        // Corrigiendo posiciones para escribirlas en el archivo
        corriger(r_x, sigma, filas);
        corriger(r_y, sigma, filas);

        //imprimo posiciones nuevas
         for(k=0; k<filas; k++)
        {
            fprintf(archivo, "%.10Lf", r_x[k]);
            fprintf(archivo, ", %.10Lf\n", r_y[k]);
        }
        fprintf(archivo, "\n");

        // Calculando energía y momento angular total
        double V = potencial(m, r_x, r_y, a_x, a_y, filas);
        double cin = cinetica(m, v_x, v_y, filas);
        double energiatotal =  cin + V;
        double momentoang = mangular(m, r_x, r_y,  v_x, v_y, filas);

        // Escribiendo en los archivos
        fprintf(archivo_, "%.10f\n", energiatotal);
        fprintf(archivo_cinetic, "%.10f\n", cin);
        fprintf(archivo_potenciale, "%.10f\n", V);
        fprintf(archivo__, "%.10f\n", momentoang);
    }

    // Liberando memoria asignada para los vectores
    free(m);
    free(r_x);
    free(r_y);
    free(v_x);
    free(v_y);
    free(w_x);
    free(w_y);
    free(a_x);
    free(a_y);

    // Cerrando archivos
    fclose(archivo);
    fclose(archivo_);
    fclose(archivo__);
    fclose(archivo_cinetic);
    fclose(archivo_potenciale);

    return 0;
}

// Función que rescala las masas
void rescm(double *masa, int n)
{
    for (int i = 0; i < n; i++)
    {
        masa[i] = masa[i] / Argonmass;
    }
}

// Función que rescata los vectores de posición
void rescr(double *pos, int sig ,int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] / sig;
    }
}

// Función que devuelve la escala original de r
void corriger(double *pos, int sig, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] * sig;
    }
}

// Función para calcular aceleración
void aceleracion(double *aceleracion, double *posx, double *posy, int L, int n)
{
    double distx, disty, distancia;
    for (int i = 0; i < n; i++)
    {
        aceleracion[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                distx = posx[i] - posx[j];
		        disty = posy[i] - posy[j];
		        distx = distx - round(distx/(2*L))*(2.*L);
		        disty = disty - round(disty/(2*L))*(2.*L);
                distancia = sqrt(distx * distx + disty * disty);
                aceleracion[i] += 24 * (2.0 / pow(distancia, 13) - 1.0 / pow(distancia, 7));
            }
        }
    }
}

// Funciones auxiliares para calcular velocidad
void calculow(double *aux, double *vel, double *acel, double h, int n)
{
    for (int i = 1; i < n; i++)
    {
        aux[i] = vel[i] + h / 2 * acel[i];
    }
}

void calculov(double *velocidad, double *aux, double *acel, double h, int n)
{
    for (int i = 1; i < n; i++)
    {
        velocidad[i] = aux[i] + h * acel[i] / 2;
    }
}

// Función para calcular la nueva posición
void calculopos(double *pos, double *vel, double *acel, double h, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] + h * vel[i] + (h * h / 2) * acel[i];
    }
}

// Función para calcular el momento angular
double mangular(double *masa, double *posx, double *posy, double *velx, double *vely, int n)
{
    double L = 0;
    for (int i = 1; i < n; i++)
    {
        L += masa[i] * sqrt(pow(posx[i], 2) + pow(posy[i], 2)) * sqrt(pow(velx[i], 2) + pow(vely[i], 2));
    }
    return L;
}

// Función para calcular la energía cinética total
double cinetica(double *masa, double *velx, double *vely, int n)
{
    double kinetic = 0;
    for (int i = 1; i < n; i++)
    {
        kinetic += 0.5 * masa[i] * (velx[i] * velx[i] + vely[i] * vely[i]);
    }
    return kinetic;
}

// Función para calcular la energía potencial total
double potencial(double *masa, double *posx, double *posy, double *acelx, double *acely, int n)
{
    double pot = 0;
    for (int i = 1; i < n; i++)
    {
        pot += -masa[i] * sqrt(pow(posx[i], 2) + pow(posy[i], 2)) * sqrt(pow(acelx[i], 2) + pow(acely[i], 2));
    }
    return pot;
}

void generar_posiciones(double *pos_x, double *pos_y, int dimension, int tamano_red)
 {
    // Número máximo de celdas en la red cuadrada
    int num_celdas = tamano_red * tamano_red;

    // Verificar si hay suficientes celdas para todas las partículas
    if (num_celdas < dimension) {
        printf("Error: No hay suficientes celdas en la red para almacenar todas las partículas.\n");
        return;
    }

    // Asignar posiciones aleatorias
    int ocupado[num_celdas]; // Array para verificar si una celda está ocupada
    for (int i = 0; i < num_celdas; i++) {
        ocupado[i] = 0; // Inicializar todas las celdas como no ocupadas
    }

    // Generar posiciones aleatorias para cada partícula
    for (int i = 0; i < dimension; i++) {
        int celda;

        // Buscar una celda no ocupada
        do {
            celda = rand() % num_celdas;
        } while (ocupado[celda]);

        // Marcar la celda como ocupada
        ocupado[celda] = 1;

        // Calcular las coordenadas (x, y) de la celda
        int fila = celda / tamano_red;
        int columna = celda % tamano_red;

        // Asignar la posición de la partícula
        pos_x[i] = (double)columna + ((double)rand() / RAND_MAX); // Posición x aleatoria entre columna y columna+1
        pos_y[i] = (double)fila + ((double)rand() / RAND_MAX);    // Posición y aleatoria entre fila y fila+1
    }
}

void generar_velocidades(double *velx, double *vely, int dimension)
{
    int valor_aleatorio, valor_aleatorio2, i;

    for(i=0; i<dimension; i++){
    valor_aleatorio =   rand() % 2; // Genera 0 o 1
    valor_aleatorio2=  rand() % 2;
    velx[i] = valor_aleatorio * 2 - 1; // Convierte 0 a -1 y 1 a +1
    vely[i]= valor_aleatorio2*2 -1;
    }
}

// Función para imprimir las coordenadas (x, y) como una matriz
void imprimirCoordenadas(double x[], double y[], int n)
 {
    // Verificar si los vectores tienen la misma longitud
    if (n <= 0) {
        printf("Los vectores están vacíos.\n");
        return;
    }

    // Imprimir las coordenadas como una matriz
    printf("Coordenadas:\n");
    for (int i = 0; i < n; i++) {
        printf("(%.2f, %.2f)\n", x[i], y[i]);
    }
}

