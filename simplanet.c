//simplanetdef

//los indices del programa empiezan en uno pues la posicion 0 de los vectores esta reservada a el "sol"
//todas las funciones tienen como primer parametro el output y el resto el imput

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
#define JUP (long double)(1.898 * pow(10, 27)) // Masa de jupiter en kg
//chat gpt
#define MARS (long double)(6.417 * pow(10, 23)) // Masa de Marte en kg
#define SATURN (long double)(5.683 * pow(10, 26)) // Masa de Saturno en kg
#define URANUS (long double)(8.681 * pow(10, 25)) // Masa de Urano en kg
#define NEPTUNE (long double)(1.024 * pow(10, 26)) // Masa de Neptuno en kg
#define PLUTO (long double)(1.309 * pow(10, 22)) // Masa de Plutón en kg

#define DMER (long double)(57.9 * pow(10, 9)) // Distancia de mercurio al sol en metros
#define DVEN (long double)(108.2 * pow(10, 9)) // Distancia de venus al sol en metros
#define DTIE (long double)(149.6 * pow(10, 9)) // Distancia de la tierra al sol en metros
#define DMAR (long double)(227.9 * pow(10, 9)) // Distancia de marte al sol en metros^
#define DJUP (long double)(7.78 * pow(10, 11)) // Distancia de jupiter al sol en metros
//chatgpt
#define DMARS (long double)(227.9 * pow(10, 9)) // Distancia de Marte al Sol en metros
#define DSATURN (long double)(1.429 * pow(10, 12)) // Distancia de Saturno al Sol en metros
#define DURANUS (long double)(2.871 * pow(10, 12)) // Distancia de Urano al Sol en metros
#define DNEPTUNE (long double)(4.495 * pow(10, 12)) // Distancia de Neptuno al Sol en metros
#define DPLUTO (long double)(5.906 * pow(10, 12)) // Distancia de Plutón al Sol en metros


#define VMER (long double)(47.87 * pow(10, 3)) // Velocidad de mercurio en m/s
#define VVEN (long double)(35.02 * pow(10, 3)) // Velocidad de venus en m/s
#define VTIE (long double)(29.78 * pow(10, 3)) // Velocidad de la tierra en m/s
#define VMAR (long double)(24.077 * pow(10, 3)) // Velocidad de marte en m/s
#define VJUP (long double)(13.07 * pow(10, 3)) // Velocidad de jupiter en m/s

//chat gpt
#define VMARS (long double)(24.077 * pow(10, 3)) // Velocidad de Marte en m/s
#define VSATURN (long double)(9.69 * pow(10, 3)) // Velocidad de Saturno en m/s
#define VURANUS (long double)(6.8 * pow(10, 3)) // Velocidad de Urano en m/s
#define VNEPTUNE (long double)(5.43 * pow(10, 3)) // Velocidad de Neptuno en m/s
#define VPLUTO (long double)(6.1 * pow(10, 3)) // Velocidad de Plutón en m/s



//Declaracion de funciones
void rescm(long double *masa, int n);
void rescr(long double *pos, int n);
void reesc_v(long double *vel, int n);
void aceleracion(long double *aceleracion, long double *masa, long double *posx, long double *posy, int n);
void calculopos(long double *pos, long double *vel, long double *acel, double h, int n);
void calculov(long double *velocidad, long double *aux, long double *acel, double h, int n);
void calculow(long double *aux, long double *vel, long double *acel, double h, int n);
void energia(long double *m, long double *r_x, long double *r_y, long double *v_x, long double *v_y, int n, long double *E);
void mangular(long double *r_x, long double *r_y, long double *v_x, long double *v_y, long double *m, int n, long double *L);




int main(void) 
{
  //abrimos fichero posiciones
 FILE *archivo;
    char simplanet[] = "simplanet.txt";
    archivo = fopen(simplanet, "w");

     //abrimos fichero energia y p
 FILE *archivo_;
    char conservacion[] = "conservacion.txt";
    archivo_ = fopen(conservacion, "w");


  //definicion de variables
    int i, step, filas, j, k;
    long double *m, *r_x, *r_y, *v_x, *v_y,*w_x, *w_y ,*a_x, *a_y, *E,*L;
    double h, t;

 //numero de planetas filas=N+1 donde N es nº planetas ojo cuenta con el sol
    filas=10;
   
    //definimos
    h=0.01;
    t=10000;
    step=t/h;

    // Asignar memoria para el vector m
    m = (long double *) malloc(filas * sizeof(long double));
 
    //Asignar memoria a vectores de posicion
    r_x = (long double *) malloc(filas * sizeof(long double));
    r_y = (long double *) malloc(filas * sizeof(long double));

    //Asignar memoria a vectores velocidad 
    v_x = (long double *) malloc(filas * sizeof(long double));
    v_y = (long double *) malloc(filas * sizeof(long double));

    //Asignar memoria a vector wx wy auxiliares
    w_x = (long double *) malloc(filas * sizeof(long double));
    w_y = (long double *) malloc(filas * sizeof(long double));


     //Asignar memoria a vectores aceleración
    a_x = (long double *) malloc(filas * sizeof(long double));
    a_y = (long double *) malloc(filas * sizeof(long double));

    //Asignar memoria vector energia y momento angular
     E = (long double *) malloc(filas * sizeof(long double));
     L = (long double *) malloc(filas * sizeof(long double));

     //mio

    // Inicializar el vector m con algunas masas
  //  m[0]=  SOL;
  //  m[1] = VEN;
 //   m[2] = TIE;
  //  m[3] = JUP;

    //inicializar vector r con los radios de los planetas
  //  r_x[0]=0;
   // r_x[1]=DVEN;
   // r_x[2]=DTIE;
   // r_x[3]=DJUP;

    //inicializar componente y de la velocidad a la velocidad orbital de los planetas
   // v_y[0]=0;
   // v_y[1]=VVEN;
  //  v_y[2]=VTIE;
   // v_y[3]=VJUP;

    //chat gpt
    // Inicialización de los vectores existentes con las distancias, masas y velocidades
// r_x: Distancias de los planetas al Sol en metros
r_x[0] = 0;
r_x[1] = DMER;
r_x[2] = DVEN;
r_x[3] = DTIE;
r_x[4] = DMAR;
r_x[5] = DJUP;
r_x[6] = DSATURN;
r_x[7] = DURANUS;
r_x[8] = DNEPTUNE;
r_x[9] = DPLUTO;

// m: Masas de los planetas en kilogramos
m[0] = SOL;
m[1] = MER;
m[2] = VEN;
m[3] = TIE;
m[4] = MAR;
m[5] = JUP;
m[6] = SATURN;
m[7] = URANUS;
m[8] = NEPTUNE;
m[9] = PLUTO;

// v_y: Velocidades en la componente y de los planetas en m/s
v_y[0] = 0;
v_y[1] = VMER;
v_y[2] = VVEN;
v_y[3] = VTIE;
v_y[4] = VMAR;
v_y[5] = VJUP;
v_y[6] = VSATURN;
v_y[7] = VURANUS;
v_y[8] = VNEPTUNE;
v_y[9] = VPLUTO;

reesc_v(v_y, filas);

    //inicializar la componente "y" y el vector "a" y "v" a 0
    for(i=0; i<filas; i++)
    {
        r_y[i]=0;
        v_x[i]=0;
       // v_y[i]=0;
        a_x[i]=0;
        a_y[i]=0; 
        E[i]=0;
        L[i]=0;
        w_x[i]=0;
        w_y[i]=0;
     }


    // Llamar a la función para reescalar el vector m rx y ry
    rescm(m, filas);
    rescr(r_x, filas);
    rescr(r_y, filas);
    
    //llamar a la aceleracion para la primera iteracion
    //calculo la nueva aceleración para x
    aceleracion(a_x, m, r_x, r_y, filas);
    //calculo la nueva aceleración para y
    aceleracion(a_y, m, r_y, r_x, filas);

    //llamar a verlet
    for(int a=0; a<step; a++)
    {

       //calculo nueva posicion x e y
        calculopos(r_x, v_x, a_x, h, filas);
        calculopos(r_y, v_y, a_y, h, filas);


        //imprimo posiciones nuevas
         for(k=0; k<filas; k++)
        {
            fprintf(archivo, "%.10Lf", r_x[k]);
            fprintf(archivo, ", %.10Lf\n", r_y[k]);
        }
        fprintf(archivo, "\n");

        //calculo w
        calculow(w_x,v_x,a_x, h, filas);
        calculow(w_y,v_y,a_y, h, filas);



        //calculo la nueva aceleración para x
        aceleracion(a_x, m, r_x, r_y, filas);
        //calculo la nueva aceleración para y
        aceleracion(a_y, m, r_y, r_x, filas);


        //calculo la nueva velocidad
        calculov(v_x , w_x, a_x, h, filas);
        calculov(v_y , w_y, a_y, h, filas);

    }

       // temp=h+temp;


        //inicializo energias
       // inic0(ec);
        //inic0(ep);

        //Sacamos el output al archivo de energia
       

     //   for(k=0; k<n; k++)
      //  {
           //fprintf(archivo, "%.10Lf %.10Lf,\n", r_x[k], r_y[k]);
        //    fprintf(archivo_, "%.10Lf", E[k]);
          //  fprintf(archivo_, ", %.10Lf\n", L[k]);
       // }
       
      //  fprintf(archivo_, "\n");
     
    


    // Liberar la memoria asignada para los vectores
    free(m);
    free(r_x);
    free(r_y);
    free(v_x);
    free(v_y);
    free(a_x);
    free(a_y);
    free(E);
    free(L);

    //cerramos fichero
 fclose(archivo);
 fclose(archivo_);

    return 0;
}



//funcion reescalar masa
void rescm(long double *masa, int n)
 {
    int i;
    for (i = 0; i < n; i++) {
        masa[i] = masa[i] / SOL;
    }
}

//funcion reescalar vector posicion input vector filas output vector reescalado
void rescr(long double *pos, int n)
{
int i,j;
    for (i = 0; i < n; i++) {
        pos[i]=pos[i]/c;
    }

}

void reesc_v(long double *vel, int n)
{
    int i;

    for (i=0; i<n; i++)
    {
        vel[i]=vel[i]*pow(((c)/(G*SOL)),0.5);
      
    }
    return;
}

//funcion aceleracion imput; masa, posicion de la aceleracion que se va a calcular, la otra posicon, filas. output; aceleracion que se va a calcular
void aceleracion(long double *aceleracion, long double *masa, long double *posx, long double *posy, int n)
{
int i,j;
long double distancia, coordx, coordy;

    for (i = 0; i < n; i++)
    {
        aceleracion[i]=0;
    }
    
    for (i = 1; i < n; i++){
    for (j = 0; j < n; j++){
    
        if(i!=j)
         {
        //calculamos distancia entre planeta i y planeta j
        coordx=posx[i] - posx[j];
        coordy=posy[i] - posy[j];
        distancia=sqrt(coordx*coordx + coordy*coordy);

        //calculamos aceleracion
        aceleracion[i]+= -masa[j]*(posx[i]-posx[j])/fabs(pow(distancia, 3));
         }                    
    }
    }

}

//funcion auxiliar para actulizar velocidad. Recibe la velocidad y aceleracion que se van a auxiliar ademas de h y devuelve el auxiliar
void calculow(long double *aux, long double *vel, long double *acel, double h, int n)
{ 
    for(int i=1; i<n; i++)
    {
         aux[i]=vel[i] + h/2*acel[i];
    }
} 

//funcion que calcula la nueva velocidad, recibe el valor auxiliar calculado antes, la aceleracion y h
void calculov(long double *velocidad, long double *aux, long double *acel, double h, int n)
{
    for(int i=1; i<n; i++)
    {
      velocidad[i]=aux[i] + h*acel[i]/2;
    }
}

//funcion que actualiza la posicion, recibe vector posicion que se va a actualizar, vector velocidad, aceleracion de esa posicion y h
void calculopos(long double *pos, long double *vel, long double *acel, double h, int n)
{

    for(int i=0; i<n; i++)
    {
        pos[i]= pos[i] + h*vel[i] + (h*h/2)*acel[i];
    }
}

void energia(long double *m,long double *r_x,long double *r_y,long double *v_x,long double *v_y, int n, long double *E)
{

     for (int i = 1; i < n; i++) {
        // calculo de la energía cinética
        long double vtotal = v_x[i] * v_x[i] + v_y[i] * v_y[i];
        E[i] += 0.5 * m[i] * vtotal;

        // calculo de la energía potencial
        long double distancia = sqrt(r_x[i] * r_x[i] + r_y[i] * r_y[i]);
        E[i] -= G * SOL * m[i] / distancia;
    }
}

void mangular(long double *r_x, long double *r_y, long double *v_x, long double *v_y , long double *m, int n, long double *L)
{
for (int i = 1; i < n; i++)
{
      L[i]+=m[i]*(r_x[i]*v_y[i]-r_y[i]*v_x[i]);
}

}

