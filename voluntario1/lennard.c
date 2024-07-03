#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define Argonmass (6.6335209 * pow(10, -26)) // Masa en kg
#define Armstrong (1 * pow(10, -10))          // Armstrong en metros
#define PI (3.1415926535897932384)

// Declaración de funciones
void aceleracion(double *aceleracion, double *posx, double *posy, int L, int n);
void calculopos(double *pos, double *vel, double *acel, double h, int n);
void calculov(double *velocidad, double *aux, double *acel, double h, int n);
void calculow(double *aux, double *vel, double *acel, double h, int n);
double cinetica(double *masa, double *velx, double *vely, int n);
double potencial(double *masa, double *posx, double *posy, double *acelx, double *acely, int n);
void generar_posiciones(double *posx, double *posy, int dimension, int longitud);
void generar_velocidades(double *velx, double *vely, int dimension);
void imprimirCoordenadas(double *x, double *y, int n);
void contorno(double *pos, int tampart, int longitud);
double temperatura(double *vx, double *vy, int final, int inicial, int dim);
void generate_array(double *posx, double *posy, int n, int L);

int main(void)
{
    // Abriendo ficheros
    FILE *archivo;
    char simplanet[] = "lennard.txt";
    archivo = fopen(simplanet, "w");

    FILE *velocity;
    char veloz[] = "velocity.txt";
    velocity = fopen(veloz, "w");

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
    char tempo[] = "tiempo.txt";
    archivo__ = fopen(tempo, "w");

    // Definición de variables
    int i, filas, j, k, N, L, sigma, epsilon;
    double *m, *r_x, *r_y, *v_x, *w_x, *w_y, *v_y, *a_x, *a_y, t, h, velocidades, T;

    t = 0;
    h = 0.002; // paso
    sigma = 1.;
    epsilon = 1.;

    L = 10; // tamaño de la caja
    N = 15; // numero de particulas
    filas = N;

    // Asignando memoria para los vectores
    m = (double *)malloc((filas+1) * sizeof(double));
    r_x = (double *)malloc((filas+1) * sizeof(double));
    r_y = (double *)malloc((filas+1) * sizeof(double));
    v_x = (double *)malloc((filas+1) * sizeof(double));
    v_y = (double *)malloc((filas+1) * sizeof(double));
    w_x = (double *)malloc((filas+1) * sizeof(double));
    w_y = (double *)malloc((filas+1) * sizeof(double));
    a_x = (double *)malloc((filas+1) * sizeof(double));
    a_y = (double *)malloc((filas+1) * sizeof(double));

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
    // generar_posiciones(r_x, r_y, N, L);
    generate_array(r_x, r_y, N, L);

    //generamos velocidades
    generar_velocidades(v_x, v_y, filas);

    //mostramos vectores posiciones por pantalla a ver si esta bien
    imprimirCoordenadas(r_x, r_y, filas);

    // Calculando aceleraciones para la primera iteración
    aceleracion(a_x, r_x, r_y, L, filas);
    aceleracion(a_y, r_y, r_x, L, filas);

    /* ALGORITMO DE VERLET */
    for (int paso = 0; paso < 100000; paso++)
    {

        // Calculando nuevas posiciones
        calculopos(r_x, v_x, a_x, h, filas);
        calculopos(r_y, v_y, a_y, h, filas);

        //aplicamos condiciones de contorno
        contorno(r_x, N, L);
        contorno(r_y, N, L);

        // Calculando w
        calculow(w_x, v_x, a_x, h, filas);
        calculow(w_y, v_y, a_y, h, filas);

        // Calculando nuevas aceleraciones
        aceleracion(a_x, r_x, r_y, L, filas);
        aceleracion(a_y, r_y, r_x, L, filas);

        // Calculando nuevas velocidades
        calculov(v_x, w_x, a_x, h, filas);
        calculov(v_y, w_y, a_y, h, filas);

        if(paso%1000==0){
        //imprimo posiciones nuevas
        for (k = 0; k < filas; k++)
        {
            fprintf(archivo, "%.10f", r_x[k]);
            fprintf(archivo, ", %.10f\n", r_y[k]);
        }
        fprintf(archivo, "\n");
        }

        // Calculando energía y momento angular total
        double V = potencial(m, r_x, r_y, a_x, a_y, filas);
        double cin = cinetica(m, v_x, v_y, filas);
        double energiatotal = cin + V;

    
        //calculamos la temperatura
        int inicial=20;
        int final=50;
        velocidades=0.;

        if(paso>=inicial && paso<=final){
        for(int i=0; i<filas; i++)
        {
        velocidades+= v_x[i]*v_x[i] + v_y[i]*v_y[i];
        }
        }

        if(paso==final){
         T= 0.5*velocidades/(filas*1.0*(final-inicial));
        printf ("%.10f\n", T);
        }

        

    /*   //para el histograma de velocidades:
        if (paso >= paso/2 && paso % 10 == 0) { // Guardar velocidades después de alcanzar el equilibrio
            for (int i = 0; i < filas; i++) {
                fprintf(velocity, "%.5f %.5f\n", v_x[i], v_y[i]);
            }
        }
     */ 
        // Escribiendo en los archivos
        fprintf(archivo_, "%.10f\n", energiatotal);
        fprintf(archivo_cinetic, "%.10f\n", cin);
        fprintf(archivo_potenciale, "%.10f\n", V);

         fprintf(archivo__, "%i\n", paso);
        
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
    fclose(velocity);
    fclose(archivo_);
    fclose(archivo__);
    fclose(archivo_cinetic);
    fclose(archivo_potenciale);

    return 0;
}

// Función que calcula las aceleraciones
void aceleracion(double *aceleracion, double *posx,double *posy, int L, int n)
{
    long double distx, disty, distancia, aux1, aux2;
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
/*
                distx = fabs(posx[i] - posx[j]);
                disty = fabs(posy[i] - posy[j]);
                aux1= L - distx;
                aux2= L - disty;

                if(distx>aux1) distx=aux1;
                else distx=distx;

                if(disty>aux2) disty=aux2;
                else disty=disty;

                */
 /*               double dx = fabs(x2 - x1);
    double dy = fabs(y2 - y1);

    // Aplicar condiciones periódicas
    if (dx > L / 2.0) dx = L - dx;
    if (dy > L / 2.0) dy = L - dy;

    */

   /*
                distx = (posx[i] - posx[j]) -L*round((posx[i] - posx[j])/L);
                disty = (posy[i] - posy[j]) -L*round((posy[i] - posy[j])/L);

    */
    
                distx = fabs(posx[i] - posx[j]);
                disty = fabs(posy[i] - posy[j]);
               if (distx > L / 2.0) distx = L - distx;
                else if (distx < L/2.) distx = L + distx;
                if (disty > L / 2.0) disty = L - disty;
                else if (distx < L/2.) disty = L + disty;
                
               // distx = distx - round(distx / (2 * L)) * (2. * L);
               // disty = disty - round(disty / (2 * L)) * (2. * L);
                distancia = sqrt(distx * distx + disty * disty);
               //  aceleracion[i] += ((posx[i] - posx[j]) / distancia) * 24. * (2.0 / pow(distancia, 13) - 1.0 / pow(distancia, 7));
                aceleracion[i] +=  ((posx[i] - posx[j])/distancia)*(48./ pow(distancia, 13)) - (24.0 / pow(distancia, 7)) ;
            }
        }
    }
}

// Función que calcula las nuevas posiciones
void calculopos(double *pos, double *vel, double *acel, double h, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] += h * vel[i] + ((h * h * acel[i]) / 2);
    }
}

// Función que calcula las nuevas velocidades
void calculov(double *velocidad, double *aux, double *acel, double h, int n)
{
    for (int i = 0; i < n; i++)
    {
        velocidad[i] = aux[i] + (h * acel[i] / 2.);
    }
}

// Función que calcula w
void calculow(double *aux, double *vel, double *acel, double h, int n)
{
    for (int i = 0; i < n; i++)
    {
        aux[i] = vel[i] + (h * acel[i] / 2.);
    }
}

// Función que calcula la energía cinética
double cinetica(double *masa, double *velx, double *vely, int n)
{
    int i;
    double ecinetica = 0.0;

     for ( i = 0; i < n; i++)
    {
         masa[i]=1.;
    }

    for ( i = 0; i < n; i++)
    {
        ecinetica += masa[i] * (velx[i] * velx[i] + vely[i] * vely[i]);
    }

    ecinetica = 0.5 * ecinetica;

    return ecinetica;
}

// Función que calcula la energía potencial
double potencial(double *masa, double *posx, double *posy, double *acelx, double *acely, int n)
{
    double epotencial = 0.0;
    double dx, dy, r, rc, sc;
    rc = pow(2, 1. / 6.);
    sc = 1;

    for (int i = 0; i < n - 1; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            dx = (posx[j] - posx[i]);
            dy = (posy[j] - posy[i]);
            r = sqrt(dx * dx + dy * dy);

            if (r <= rc)
            {
                epotencial += 4 * (pow(r, -12) - pow(r, -6)) + 1;
            }
        }
    }

    return epotencial;
}

// Función que genera posiciones aleatorias
void generar_posiciones(double *posx, double *posy, int dimension, int longitud)
{
    for (int i = 0; i < dimension; i++)
    {
        posx[i] = ((double)rand() / (double)(RAND_MAX)) * longitud;
        posy[i] = ((double)rand() / (double)(RAND_MAX)) * longitud;
    }
}

// Función que genera velocidades aleatorias
void generar_velocidades(double *velx, double *vely, int dimension)
{
    int valor_aleatorio, valor_aleatorio2, i;
    double theta, r;
    r = 1.;

    for (i = 0; i < dimension; i++)
    {
        theta = 2. * PI * rand() / RAND_MAX;
        velx[i] = r * cos(theta);
        vely[i] = r * sin(theta);
    }
}

// Función que imprime las coordenadas
void imprimirCoordenadas(double *x, double *y, int n)
{
    for (int i = 0; i < n; i++)
    {
        printf("x: %.10f, y: %.10f\n", x[i], y[i]);
    }
}


void contorno(double *pos, int tampart, int longitud) 
{
    int i;
    long double k;

    for (i = 0; i < tampart; i++) {
     //  k=fmod(fabs(pos[i]), longitud) +0.0;
        if (pos[i] > longitud) 
        {  // if(k==0) pos[i] = (long double);
            pos[i] =  fmod(pos[i], longitud);

        } else if (pos[i] < 0.) {
            pos[i] =  longitud - fmod(fabs(pos[i]), longitud);
        }
        else pos[i]=pos[i];
    }
}

// Función que calcula la temperatura
double temperatura(double *vx, double *vy, int final, int inicial, int dim)
{
    double temp = 0.0;
    int count = 0;

    for (int i = inicial; i < final; i++)
    {
        temp += vx[i] * vx[i] + vy[i] * vy[i];
        count++;
    }

    temp = temp / (2.0 * count);

    return temp;
}

// Genera un array con posiciones aleatorias para x e y
void generate_array( double *posx, double *posy, int n, int L) 
{
 //   srand(time(NULL));  // Inicializa el generador de números aleatorios

    for (int i = 0; i < n; i++) {
        int value1 = (i % (L / 2)) * 2 + 1;
        int value2 = (i % 4) * 2 + 1;
        
        // Genera dos números aleatorios entre 0 y 1
        float rand1 = (float)rand() / RAND_MAX;
        float rand2 = (float)rand() / RAND_MAX;

        // Asigna los valores calculados más los aleatorios a los punteros posx y posy
        posx[i] = (double) 1.0*(value1 + rand1);
        posy[i] = ( double) 1.0*(value2 + rand2);
    }
}