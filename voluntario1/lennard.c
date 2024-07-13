#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define Argonmass (6.6335209 * pow(10, -26)) // Masa en kg
#define Armstrong (1 * pow(10, -10))          // Armstrong en metros
#define PI (3.1415926535897932384)

// Declaración de funciones
void aceleracion(double *aceleracionx, double *aceleraciony, double *posx, double *posy, int L, int n);
void calculopos(double *pos, double *vel, double *acel, double h, int n);
void calculov(double *velocidad, double *aux, double *acel, double h, int n);
void calculow(double *aux, double *vel, double *acel, double h, int n);
double cinetica(double *masa, double *velx, double *vely, int n);
double potencial(double *masa, double *posx, double *posy, int n, int L);
void generar_posiciones(double *posx, double *posy, int dimension, int longitud);
void generar_velocidades(double *velx, double *vely, int dimension);
void imprimirCoordenadas(double *x, double *y, int n);
void contorno(double pos[] , double velocidad[], double *momento, int tampart, int longitud);
double temperatura(double *vx, double *vy, int final, int inicial, int dim);
void generate_array(double *posx, double *posy, int n, int L);
void posicioncuad(double *posx, double *posy, int n, int L);
double aleatorio();
void calentarrapido(double *vx, double *vy, int N, double aumento);
double distanciapart(double posx1, double posx2, double posy1,double posy2, int N, int L);

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

    FILE *archivuo;
    char fluctu[] = "fluctuaciones.txt";
    archivuo = fopen(fluctu, "w");

     FILE *temp;
    char tempera[] = "temperatura.txt";
    temp = fopen(tempera, "w");

    FILE *correlacion;
    char corre[] = "correlacion.txt";
    correlacion = fopen(corre, "w");


    // Definición de variables
    int i, filas, j, k, N, L, sigma, epsilon, numiter, inicio, fin, veces, Lindice;
    double *m, *r_x, *r_y, *v_x, *w_x, *w_y, *v_y, *a_x, *a_y, *momento, t, h, velocidades, T, suma, vesp,momentillo, P, aumentoT, *r_xo, *r_yo, flux, Ttiempo;
    double xfijo, yfijo, xvariab, yvariab, correcaminos;

    t = 0;
    h = 0.002; // paso
    sigma = 1.;
    epsilon = 1.;

    L = 4; // tamaño de la caja
    N = 16; // numero de particulas
    filas = N;
    numiter=2.*50000; //numero de iteraciones del gas
    inicio=20; //inicio para histograma velocidades
    fin=50; //fin para histograma velocidades
    momentillo=0.0; // inicializamos el momento
    Lindice=0;
    aumentoT=1.1;
    flux=0.;

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
    momento = (double *)malloc((filas+1) * sizeof(double));
    r_xo = (double *)malloc((filas+1) * sizeof(double));
    r_yo = (double *)malloc((filas+1) * sizeof(double));

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

    veces=0;
    suma=0.0;

    //Asignamos posiciones a las particulas
    // generar_posiciones(r_x, r_y, N, L);
    //generate_array(r_x, r_y, N, L);
    posicioncuad(r_x, r_y, N, L);

    //apartado de calentar
    for(i=0; i<N; i++)
    {
        r_xo[i]=r_x[i];
        r_yo[i]=r_y[i];
    }

    //generamos velocidades
    generar_velocidades(v_x, v_y, filas);

    //mostramos vectores posiciones por pantalla a ver si esta bien
    imprimirCoordenadas(r_x, r_y, filas);

    // Calculando aceleraciones para la primera iteración
    aceleracion(a_x, a_y, r_x, r_y, L, filas);

    /* ALGORITMO DE VERLET */
    for (int paso = 0; paso < numiter; paso++)
    {

        // Calculando nuevas posiciones
        calculopos(r_x, v_x, a_x, h, filas);
        calculopos(r_y, v_y, a_y, h, filas);

        //aplicamos condiciones de contorno
        contorno(r_x, v_x, momento, N, L);
        contorno(r_y, v_y, momento, N, L);

        // Calculando w
        calculow(w_x, v_x, a_x, h, filas);
        calculow(w_y, v_y, a_y, h, filas);

        // Calculando nuevas aceleraciones
        aceleracion(a_x, a_y, r_x, r_y, L, filas);

        // Calculando nuevas velocidades
        calculov(v_x, w_x, a_x, h, filas);
        calculov(v_y, w_y, a_y, h, filas);

       // if(paso*h==20 || paso*h==30 || paso*h==35 || paso*h==45){
       if(paso*h==180){
        calentarrapido(v_x, v_y, N, aumentoT);
        }

        Ttiempo=0.;
        for(int g=0; g<N; g++){
            Ttiempo+=v_x[g]*v_x[g] + v_y[g]*v_y[g];
        }
        Ttiempo=Ttiempo/(N);

        fprintf(temp, "%.10f", Ttiempo);
         fprintf(temp, "\n");

        flux=0.;
        for(int s=0; s<N; s++)
        {
            flux+= (r_x[s]+ r_y[s]-r_xo[s]-r_yo[s])*(r_x[s]+ r_y[s]-r_xo[s]-r_yo[s]);
        }
        flux=flux/N;
        
      // flux+=(r_x[0]+ r_y[0]-r_xo[0]-r_yo[0])*(r_x[0]+ r_y[0]-r_xo[0]-r_yo[0]);

         fprintf(archivuo, "%.10f", flux);
         fprintf(archivuo, "\n");

    //almacenamos la posicion de la particula cada x iteraciones porque si no se ve muy lento
          if(paso%1000==0){
       // if(paso%100==0){
        //imprimo posiciones nuevas
        for (k = 0; k < filas; k++)
        {
            fprintf(archivo, "%.10f", r_x[k]);
            fprintf(archivo, ", %.10f\n", r_y[k]);
        }
        fprintf(archivo, "\n");
        }

        correcaminos=0.;
          // El de correlación de pares
         /* for(i=0; i<N; i++){
               if(i!=0) {
                    xfijo=r_x[0];
                    yfijo=r_y[0];
                    xvariab=r_x[i];
                    yvariab=r_y[i];
                  correcaminos=  distanciapart(xfijo, xvariab, yfijo, yvariab, N, L);

                     fprintf(correlacion, "%.10f\n", correcaminos);
                 
                }
                
            }
             */
            
                    xfijo=r_x[0];
                    yfijo=r_y[0];
                    xvariab=r_x[5];
                    yvariab=r_y[5];
                  correcaminos=  distanciapart(xfijo, xvariab, yfijo, yvariab, N, L)*distanciapart(xfijo, xvariab, yfijo, yvariab, N, L);

            fprintf(correlacion, "%.10f\n", correcaminos);


        // Calculando energía y momento angular total
        double V = potencial(m, r_x, r_y, N, L);
        double cin = cinetica(m, v_x, v_y, filas);
        double energiatotal = cin + V;

    
       //para el histograma de velocidades:
        if (paso*h >=inicio && (paso*h) <=fin) { // Guardar velocidades después de alcanzar el equilibrio
            for (int i = 0; i < filas; i++) {
     //           fprintf(velocity, "%.5f %.5f %.5f\n", v_x[i], v_y[i], v_x[i]*v_x[i] + v_y[i]*v_y[i]); //sin modulo de velocidades
                 fprintf(velocity, "%.5f  %.5f %.5f %.5f\n", v_x[i], v_y[i], sqrt(v_x[i]*v_x[i] + v_y[i]*v_y[i]), fabs(v_x[i])); //con modulo de velocidades

                 suma+=v_x[i]*v_x[i] + v_y[i]*v_y[i];
                 veces++;
            }
            momentillo+=momento[0]/(paso*h*L*L); 
            Lindice++;
        }
      
        // Escribiendo en los archivos
        fprintf(archivo_, "%.10f\n", energiatotal);
        fprintf(archivo_cinetic, "%.10f\n", cin);
        fprintf(archivo_potenciale, "%.10f\n", V);

         fprintf(archivo__, "%f\n", paso*h);
        
    }

    T=suma/(2.0*veces);
    P=momentillo/Lindice;

    vesp=sqrt(T);

    printf("La temperatura es: ");
    printf("%.10f\n", T);
    printf("La velocidad esperada es: ");
    printf("%.10f\n", vesp);
    printf("La presion es: ");
    printf("%.10f\n", P);

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
    free(momento);
    free(r_xo);
    free(r_yo);

    // Cerrando archivos
    fclose(archivo);
    fclose(velocity);
    fclose(archivo_);
    fclose(archivo__);
    fclose(archivo_cinetic);
    fclose(archivo_potenciale);
    fclose(archivuo);
    fclose(temp);
    fclose(correlacion);

    return 0;
}

// Función que calcula las aceleraciones
void aceleracion(double *aceleracionx, double *aceleraciony, double *posx, double *posy, int L, int n)
{
    double distx, disty, distancia, aux1, aux2, rdif[2];
    for (int i = 0; i < n; i++)
    {
        aceleracionx[i] = 0;
        aceleraciony[i] = 0;
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                // Obtengo el cuadrado determinado por la retícula de la partícula 2 en cuyo interior se encuentra la partícula 1
            double incx=posx[i]-posx[j];
            double incxraro;

            if(posx[i]<=posx[j]) incxraro=posx[i]-posx[j]+L;
            else incxraro=posx[i]-posx[j]-L;

            double incy=posy[i]-posy[j];
            double incyraro;

            if(posy[i]<=posy[j]) incyraro=posy[i]-posy[j]+L;
            else incyraro=posy[i]-posy[j]-L;

            // Calculo las posibles distancias al cuadrado, la que sea mínima es la verdadera
            double d1=incx*incx+incy*incy;
            double d2=incxraro*incxraro+incy*incy;
            double d3=incx*incx+incyraro*incyraro;
            double d4=incxraro*incxraro+incyraro*incyraro;

            // Si la distancia usual es la menor, el vector es el usual
            if(d1<=d2 && d1<=d3 && d1<=d4) {
            rdif[0]=incx;
            rdif[1]=incy;
            distancia= sqrt(d1);
            }

            // Si es la segunda, pues el correspondiente, y así
            else if(d2<=d1 && d2<=d3 && d2<=d4) {
            rdif[0]=incxraro;
            rdif[1]=incy;
            distancia = sqrt(d2);
            }

            else if(d3<=d1 && d3<=d2 && d3<=d4) {
            rdif[0]=incx;
            rdif[1]=incyraro;
            distancia = sqrt(d3);
            }

            else {
            rdif[0]=incxraro;
            rdif[1]=incyraro;
            distancia= sqrt(d4);
            }

                if(distancia<3.){
              //  aceleracion[i] +=  ((posx[i] - posx[j])/distancia)*(48./ pow(distancia, 13)) - (24.0 / pow(distancia, 7)) ;
                aceleracionx[i]+=24*(2-pow(distancia,6))*rdif[0]/(pow(distancia,14));
                aceleraciony[i]+=24*(2-pow(distancia,6))*rdif[1]/(pow(distancia,14));
                }
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
double potencial(double *masa, double *posx, double *posy, int n, int L)
{
    int i,j;
    double epotencial = 0.0;
    double distx, disty, distancia;

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
// Obtengo el cuadrado determinado por la retícula de la partícula 2 en cuyo interior se encuentra la partícula 1
            double incx=posx[i]-posx[j];
            double incxraro;

            if(posx[i]<=posx[j]) incxraro=posx[i]-posx[j]+L;
            else incxraro=posx[i]-posx[j]-L;

            double incy=posy[i]-posy[j];
            double incyraro;

            if(posy[i]<=posy[j]) incyraro=posy[i]-posy[j]+L;
            else incyraro=posy[i]-posy[j]-L;

            // Calculo las posibles distancias al cuadrado, la que sea mínima es la verdadera
            double d1=incx*incx+incy*incy;
            double d2=incxraro*incxraro+incy*incy;
            double d3=incx*incx+incyraro*incyraro;
            double d4=incxraro*incxraro+incyraro*incyraro;

            // Si la distancia usual es la menor, el vector es el usual
            if(d1<=d2 && d1<=d3 && d1<=d4) {
            distancia= sqrt(d1);
            }

            // Si es la segunda, pues el correspondiente, y así
            else if(d2<=d1 && d2<=d3 && d2<=d4) {
            distancia = sqrt(d2);
            }

            else if(d3<=d1 && d3<=d2 && d3<=d4) {
            distancia = sqrt(d3);
            }

            else {
            distancia= sqrt(d4);
            }

                epotencial +=  4 * (pow((distancia), -12) - pow((distancia), -6));
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
    r = 4.;

    for (i = 0; i < dimension; i++)
    { 
        theta = 2. * PI * rand() / RAND_MAX;
         // velx[i] = (r * cos(theta));
         // vely[i] = r * sin(theta);
         //   velx[i] = fabs(r * cos(theta));
           vely[i] = 0.;
          velx[i] = 0.;
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


void contorno(double pos[] , double velocidad[], double *momento, int tampart, int longitud) 
{
    int i;
    // Itero sobre átomos y dimensiones
      // Para cada partícula
    for(int i=0; i<tampart; i++) {

            // Si ha cruzado, añade la contribución correspondiente a momentosparaP y reduce
            if(pos[i]>= longitud || pos[i]<0) {
                momento[0]+=2.*fabs(velocidad[i]);

                pos[i]+= -longitud*floor(pos[i]/longitud);
            }
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

void posicioncuad(double *posx, double *posy, int n, int L)
{
    // Estado sólido: red cuadrada y parten del reposo
        for (int i= 0; i<n; i++){
        posx[i]=0.5+i%L;
        posy[i]=0.5+(i-i%L)/(1.0*L);
        }
}

double aleatorio() 
{
    return (double)rand() / (double)RAND_MAX;
}

void calentarrapido(double *vx, double *vy, int N, double aumento)
{
    for(int i=0; i<N; i++)
    {
        vx[i]=vx[i]*aumento;
        vy[i]=vy[i]*aumento;
    }
}

double distanciapart(double posx1, double posx2, double posy1,double posy2, int N, int L)
{

    double distancia;
    double incx=posx1-posx2;
            double incxraro;

            if(posx1<=posx2) incxraro=posx1-posx2+L;
            else incxraro=posx1-posx2-L;

            double incy=posy1-posy2;
            double incyraro;

            if(posy1<=posy2) incyraro=posy1-posy2+L;
            else incyraro=posy1-posy2-L;

            // Calculo las posibles distancias al cuadrado, la que sea mínima es la verdadera
            double d1=incx*incx+incy*incy;
            double d2=incxraro*incxraro+incy*incy;
            double d3=incx*incx+incyraro*incyraro;
            double d4=incxraro*incxraro+incyraro*incyraro;

            // Si la distancia usual es la menor, el vector es el usual
            if(d1<=d2 && d1<=d3 && d1<=d4) {
            distancia= sqrt(d1);
            }

            // Si es la segunda, pues el correspondiente, y así
            else if(d2<=d1 && d2<=d3 && d2<=d4) {
            distancia = sqrt(d2);
            }

            else if(d3<=d1 && d3<=d2 && d3<=d4) {
            distancia = sqrt(d3);
            }

            else {
            distancia= sqrt(d4);
            }

            return distancia;
}