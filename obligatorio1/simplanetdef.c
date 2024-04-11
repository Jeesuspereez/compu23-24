//simplanetdef

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define c (long double)(1.496 * pow(10, 11)) // Unidad astronómica en metros
#define G (long double)(6.674 * pow(10, -11)) // Constante de gravitación universal en m^3 kg^-1 s^-2

#define SOL (long double)(1.99 * pow(10, 30)) // Masa del sol en kg
#define MER (long double)(3.301 * pow(10, 23)) // Masa de mercurio en kg
#define VEN (long double)(4.867 * pow(10, 24)) // Masa de venus en kg
#define TIE (long double)(5.972 * pow(10, 24)) // Masa de la tierra en kg
#define MAR (long double)(6.417 * pow(10, 23)) // Masa de marte en kg

#define DMER (long double)(57.9 * pow(10, 9)) // Distancia de mercurio al sol en metros
#define DVEN (long double)(108.2 * pow(10, 9)) // Distancia de venus al sol en metros
#define DTIE (long double)(149.6 * pow(10, 9)) // Distancia de la tierra al sol en metros
#define DMAR (long double)(227.9 * pow(10, 9)) // Distancia de marte al sol en metros

#define RMER (long double)(2.4397 * pow(10, 6)) // Radio de mercurio en metros
#define RVEN (long double)(6.0518 * pow(10, 6)) // Radio de venus en metros
#define RTIE (long double)(6.371 * pow(10, 6)) // Radio de la tierra en metros
#define RMAR (long double)(3.3895 * pow(10, 6)) // Radio de marte en metros

#define VMER (long double)(47.87 * pow(10, 3)) // Velocidad de mercurio en m/s
#define VVEN (long double)(35.02 * pow(10, 3)) // Velocidad de venus en m/s
#define VTIE (long double)(29.78 * pow(10, 3)) // Velocidad de la tierra en m/s
#define VMAR (long double)(24.077 * pow(10, 3)) // Velocidad de marte en m/s



//Declaracion de funciones
void rescm(long double *m, int filas);
void rescr(long double *r_x, long double *r_y, int n);
void aceleracion(long double *m, long double *r_x, long double *r_y, long double *a_x, long double *a_y, int n);
void verlet(long double *r_x, long double *r_y, long double *v_x, long double *v_y , long double *a_x, long double *a_y, long double *m, int n, long double h);

int main(void) {
  //abrimos fichero
 FILE *archivo;
    char simplanet[] = "simplanet.txt";
    archivo = fopen(simplanet, "w");

  //definicion de variables
    int i, columnas, filas, j, k;
    long double *m, *r_x, *r_y, *v_x, *v_y, *a_x, *a_y, h,t;

    //numero de planetas n=N+1 donde N es nº planetas
    int n=3;
    filas=n;

    //definimos
    h=0.01;
    t=0;

    // Asignar memoria para el vector m
    m = (long double *) malloc(filas * sizeof(long double*));
 
    //Asignar memoria a vectores de posicion
    r_x = (long double *) malloc(filas * sizeof(long double*));
    r_y = (long double *) malloc(filas * sizeof(long double*));

    //Asignar memoria a vectores velocidad 
    v_x = (long double *) malloc(filas * sizeof(long double*));
    v_y = (long double *) malloc(filas * sizeof(long double*));


     //Asignar memoria a vectores aceleración
    a_x = (long double *) malloc(filas * sizeof(long double*));
    a_y = (long double *) malloc(filas * sizeof(long double*));


    // Inicializar el vector m con algunas masas
    m[0] = MER;
    m[1] = MAR;
    m[2] = VEN;

    //inicializar vector r con los radios de los planetas
    r_x[0]=RMER;
    r_x[1]=RMAR;
    r_x[2]=RVEN;

    //inicializar componente y de la velocidad a la velocidad orbital de los planetas
    v_y[0]=VMER;
    v_y[1]=VMAR;
    v_y[2]=VVEN;


    //inicializar la componente "y" y el vector "a" y "v" a 0
    for(i=0; i<n; i++){
        r_y[i]=0;
        v_x[i]=0;
       // v_y[i]=0;
        a_x[i]=0;
        a_y[i]=0; 
    }


    // Llamar a la función para reescalar el vector m
    rescm(m, n);
    rescr(r_x, r_y, n);
    
    //llamar a la aceleracion para la primera iteracion
    aceleracion(m, r_x, r_y, a_x, a_y, n);

    //llamar a verlet
    for(t=0; t<1000; t=t+h){
        verlet(r_x, r_y, v_x, v_y, a_x, a_y, m, n ,h);

        for(k=0; k<n; k++)
        {
           //fprintf(archivo, "%.10Lf %.10Lf,\n", r_x[k], r_y[k]);

           fprintf(archivo, "%.10Lf", r_x[k]);
            fprintf(archivo, ", %.10Lf\n", r_y[k]);
        }
         fprintf(archivo, "\n");
     
    }


    // Liberar la memoria asignada para los vectores
    free(m);
    free(r_x);
    free(r_y);
    free(v_x);
    free(v_y);
    free(a_x);
    free(a_y);

    //cerramos fichero
 fclose(archivo);

    return 0;
}





//funcion reescalar masa
void rescm(long double *m, int n) {
    int i;
    for (i = 0; i < n; i++) {
        m[i] = m[i] / SOL;
    }
}

//funcion reescalar vector posicion
void rescr(long double *r_x, long double *r_y,  int n){
int i,j;
    for (i = 0; i < n; i++) {
        r_x[i]=r_x[i]/c;
        r_y[i]=r_y[i]/c;
    }

}

//funcion aceleracion
void aceleracion(long double *m, long double *r_x, long double *r_y, long double *a_x, long double *a_y, int n)
{
int i,j;
long double distancia, coordx, coordy;

    for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
    
    if(i!=j)
    {
        //calculamos distancia entre planeta i y planeta j
        coordx=r_x[i] - r_x[j];
        coordy=r_y[i] - r_y[j];
        distancia=sqrt(coordx*coordx + coordy*coordy);
        //calculamos coordenada x
        a_x[i]+=-m[j]*(r_x[i]-r_x[j])/pow(distancia, 3);

       //calculamos coordenada y
        a_y[i]+=-m[j]*(r_y[i]-r_y[j])/pow(distancia, 3);
    }

                                 }
                                 }
}

//algortimo verlet
void verlet(long double *r_x, long double *r_y, long double *v_x, long double *v_y , long double *a_x, long double *a_y, long double *m, int n, long double h)
{
int  i;

long double *w_x, *w_y;
  //Asignar memoria a vector w auxiliar
    w_x = (long double *) malloc(n * sizeof(long double*));
    w_y = (long double *) malloc(n * sizeof(long double*));

    for (i = 0; i <n; i++)
    {
        //calculo verlet r
        r_x[i]+=  h*v_x[i] + h*h*a_x[i]/2;
        r_y[i]+= h*v_y[i] + h*h*a_y[i]/2;

        //auxiliar para v
        w_x[i]=v_x[i] + h/2*a_x[i];
        w_y[i]=v_y[i] + h/2*a_y[i];

        //actualiza a
        aceleracion(m, r_x, r_y, a_x, a_y, n);

        v_x[i]=w_x[i] + h*a_x[i]/2;
        v_y[i]=w_y[i] + h*a_y[i]/2;
    }
}

void plineal()
{
    
}
