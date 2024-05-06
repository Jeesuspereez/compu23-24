#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define Argonmass (6.6335209L * pow(10, -26)) // Masa en kg
#define Armstrong (1 * pow(10, -10))          // Armstrong en metros
#define PI (3.1415926535897932384L)

// Declaración de funciones
void rescm(long double *masa, int n);
void rescr(long double *pos, int sig, int n);
void corriger(long double *pos, int sig, int n);
void aceleracion(long double *aceleracion, long double *posx, long double *posy, int L, int n);
void calculopos(long double *pos, long double *vel, long double *acel, long double h, int n);
void calculov(long double *velocidad, long double *aux, long double *acel, long double h, int n);
void calculow(long double *aux, long double *vel, long double *acel, long double h, int n);
long double cinetica(long double *masa, long double *velx, long double *vely, int n);
long double potencial(long double *masa, long double *posx, long double *posy, long double *acelx, long double *acely, int n);
long double mangular(long double *masa, long double *posx, long double *posy, long double *velx, long double *vely, int n);
void generar_posiciones(long double *posx, long double *posy, int dimension, int longitud);
void generar_velocidades(long double *velx, long double *vely, int dimension);
void imprimirCoordenadas(long double *x, long double *y, int n);
void contorno(long double *pos, int tampart, int longitud);

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
    long double *m, *r_x, *r_y, *v_x, *w_x, *w_y, *v_y, *a_x, *a_y, t, h;

    t = 0;
    h = 0.0002; // paso
    sigma = 1;
    epsilon = 1;

    L = 10; // tamaño de la caja
    N = 20; // numero de particulas
    filas = N;

    // Asignando memoria para los vectores
    m = (long double *)malloc((filas) * sizeof(long double));
    r_x = (long double *)malloc(filas * sizeof(long double));
    r_y = (long double *)malloc(filas * sizeof(long double));
    v_x = (long double *)malloc(filas * sizeof(long double));
    v_y = (long double *)malloc(filas * sizeof(long double));
    w_x = (long double *)malloc(filas * sizeof(long double));
    w_y = (long double *)malloc(filas * sizeof(long double));
    a_x = (long double *)malloc(filas * sizeof(long double));
    a_y = (long double *)malloc(filas * sizeof(long double));

    // Inicializando aceleraciones a cero
    for (i = 0; i < filas; i++)
    {
        m[i] = 0;
        r_x[i] = 0;
        r_y[i] = 0;
        v_x[i] = 0;
        v_y[i] = 0;
        w_x[i] = 0;
        w_y[i] = 0;
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

    generar_velocidades(v_x, v_y, filas);

    // Rescalando velocidades
    /* rescr(v_y, sigma, filas);
    corriget(v_y, filas);*/

    //mostramos vectores posiciones por pantalla a ver si esta bien
    imprimirCoordenadas(r_x, r_y, filas);

    // Calculando aceleraciones para la primera iteración
    aceleracion(a_x, r_x, r_y, L, filas);
    aceleracion(a_y, r_y, r_x, L, filas);

    /* ALGORITMO DE VERLET */
    for (t = 0; t < 100; t++)
    {

        // Calculando nuevas posiciones
        calculopos(r_x, v_x, a_x, t, filas);
        calculopos(r_y, v_y, a_y, t, filas);

        //aplicamos condiciones de contorno
        contorno(r_x, N, L);
        contorno(r_y, N, L);

        // Calculando w
        calculow(w_x, v_x, a_x, t, filas);
        calculow(w_y, v_y, a_y, t, filas);

        // Calculando nuevas aceleraciones
        aceleracion(a_x, r_x, r_y, L, filas);
        aceleracion(a_y, r_y, r_x, L, filas);

        // Calculando nuevas velocidades
        calculov(v_x, w_x, a_x, t, filas);
        calculov(v_y, w_y, a_y, t, filas);

        // Corrigiendo posiciones para escribirlas en el archivo
        corriger(r_x, sigma, filas);
        corriger(r_y, sigma, filas);

        //imprimo posiciones nuevas
        for (k = 0; k < filas; k++)
        {
            fprintf(archivo, "%.10Lf", r_x[k]);
            fprintf(archivo, ", %.10Lf\n", r_y[k]);
        }
        fprintf(archivo, "\n");

        // Calculando energía y momento angular total
        long double V = potencial(m, r_x, r_y, a_x, a_y, filas);
        long double cin = cinetica(m, v_x, v_y, filas);
        long double energiatotal = cin + V;
        long double momentoang = mangular(m, r_x, r_y, v_x, v_y, filas);

        // Escribiendo en los archivos
        fprintf(archivo_, "%.10Lf\n", energiatotal);
        fprintf(archivo_cinetic, "%.10Lf\n", cin);
        fprintf(archivo_potenciale, "%.10Lf\n", V);
        fprintf(archivo__, "%.10Lf\n", momentoang);
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
void rescm(long double *masa, int n)
{
    for (int i = 0; i < n; i++)
    {
        masa[i] = masa[i] / Argonmass;
    }
}

// Función que rescata los vectores de posición
void rescr(long double *pos, int sig, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] / sig;
    }
}

// Función que devuelve la escala original de r
void corriger(long double *pos, int sig, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] * sig;
    }
}

void aceleracion(long double *aceleracion, long double *posx, long double *posy, int L, int n)
{
    long double distx, disty, distancia;
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
               // distx = distx - round(distx / (2 * L)) * (2. * L);
               // disty = disty - round(disty / (2 * L)) * (2. * L);
                distancia = sqrt(distx * distx + disty * disty);
                aceleracion[i] += ((posx[i] - posx[j]) / distancia) * 24 * (2.0 / pow(distancia, 13) - 1.0 / pow(distancia, 7));
            }
        }
    }
}

// Funciones auxiliares para calcular velocidad
void calculow(long double *aux, long double *vel, long double *acel, long double h, int n)
{
    for (int i = 1; i < n; i++)
    {
        aux[i] = vel[i] + h / 2 * acel[i];
    }
}

void calculov(long double *velocidad, long double *aux, long double *acel, long double h, int n)
{
    for (int i = 1; i < n; i++)
    {
        velocidad[i] = aux[i] + h * acel[i] / 2;
    }
}

// Función para calcular la nueva posición
void calculopos(long double *pos, long double *vel, long double *acel, long double h, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] + h * vel[i] + (h * h / 2) * acel[i];
    }
}

// Función para calcular el momento angular
long double mangular(long double *masa, long double *posx, long double *posy, long double *velx, long double *vely, int n)
{
    long double L = 0;
    for (int i = 1; i < n; i++)
    {
        L += masa[i] * sqrt(pow(posx[i], 2) + pow(posy[i], 2)) * sqrt(pow(velx[i], 2) + pow(vely[i], 2));
    }
    return L;
}

// Función para calcular la energía cinética total
long double cinetica(long double *masa, long double *velx, long double *vely, int n)
{
    long double kinetic = 0;
    for (int i = 1; i < n; i++)
    {
        kinetic += 0.5 * masa[i] * (velx[i] * velx[i] + vely[i] * vely[i]);
    }
    return kinetic;
}

// Función para calcular la energía potencial total
long double potencial(long double *masa, long double *posx, long double *posy, long double *acelx, long double *acely, int n)
{
    long double pot = 0;
    for (int i = 1; i < n; i++)
    {
        pot += -masa[i] * sqrt(pow(posx[i], 2) + pow(posy[i], 2)) * sqrt(pow(acelx[i], 2) + pow(acely[i], 2));
    }
    return pot;
}

void generar_posiciones(long double *posx, long double *posy, int dimension, int longitud)
{
    int n, m, i, j;

    for (i = 0; i < dimension; i++)
    {
        do
        {
            posx[i] = ((long double)rand() / RAND_MAX) * (longitud - 2) + 1;
            posy[i] = ((long double)rand() / RAND_MAX) * (longitud - 2) + 1;

            for (j = 0; j < i; j++)
            {
                long double distancia_x = posx[i] - posx[j];
                long double distancia_y = posy[i] - posy[j];
                long double distancia_entre_particulas = sqrt(distancia_x * distancia_x + distancia_y * distancia_y);
                if (distancia_entre_particulas < 1.7)
                {
                    break;
                }
            }
        } while (j < i);
    }
}

void generar_velocidades(long double *velx, long double *vely, int dimension)
{
    int valor_aleatorio, valor_aleatorio2, i;
    long double theta, r;
    r = 1.;

    for (i = 0; i < dimension; i++)
    {
        theta = 2. * PI * rand() / RAND_MAX;
        velx[i] = r * cos(theta);
        vely[i] = r * sin(theta);
    }
}

// Función para imprimir las coordenadas (x, y) como una matriz
void imprimirCoordenadas(long double *x, long double *y, int n)
{
    // Imprimir las coordenadas como una matriz
    printf("Coordenadas:\n");
    for (int i = 0; i < n; i++)
    {
        printf("(%.2Lf, %.2Lf)\n", x[i], y[i]);
    }
}

void contorno(long double *pos, int tampart, int longitud)
{
    int i, k;

    for (i = 0; i < tampart; i++)
    {
        k = 0;
        if (pos[i] > longitud)
        {
            k = floor(pos[i] / longitud); //la funcion floor devuelve el resto real de la fraccion por ej pos=2.5L --> floor(pos/long)=.5
            pos[i] = pos[i] - k * longitud;
        }
        else if (pos[i] < 0)
        {
            k = floor(-pos[i] / longitud);
            pos[i] = longitud + pos[i] + k * longitud;
        }
        else
            pos[i] = pos[i];
    }
}
