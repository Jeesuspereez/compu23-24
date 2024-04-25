/* programa que simula la dinámica molecular de un gas con un potencial de lennard jones
-V(r)=4eps[(sigma/r)^12 - (sigma/r)^6]
-F(r)= -parcialV/r = 24eps*sigma^6(2sigma^6-r^6)/r^13 */  

//los indices del programa empiezan en uno pues la posicion 0 de los vectores esta reservada a el "sol"
//todas las funciones tienen como primer parametro el output y el resto el imput

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//libreria que optimiza bucles
#include <omp.h>

#define Argonmass (long double)(6.6335209 * pow(10, -26)) // Masa en kg
#define Armstrong (long double)(1 * pow(10, -10)) // Armstrong en metros

//Declaracion de funciones
void rescm(long double *masa, int n);
void rescr(long double *pos, int n);
void corriger(long double *pos, int n);
void reesct(long double *tiempo, int n);
void corriget(long double *tiempo, int n);
void aceleracion(long double *aceleracion, long double *masa, long double *posx, long double *posy, int n);
void calculopos(long double *pos, long double *vel, long double *acel, double h, int n);
void calculov(long double *velocidad, long double *aux, long double *acel, double h, int n);
void calculow(long double *aux, long double *vel, long double *acel, double h, int n);
long double cinetica(long double *masa, long double *velx, long double *vely, int n);
long double potencial(long double *masa, long double *posx, long double *posy, long double *acelx, long double *acely,int n);
long double mangular(long double *masa, long double *posx, long double *posy, long double *velx, long double *vely, int n);

int main(void) 
{
  //abrimos fichero posiciones
 FILE *archivo;
    char simplanet[] = "nplanet.txt";
    archivo = fopen(simplanet, "w");

     //abrimos fichero cinetica
    FILE *archivo_cinetic;
    char cinetic[] = "cinetica.txt";
    archivo_cinetic = fopen(cinetic, "w");

    //abrimos fichero potencial
    FILE *archivo_potenciale;
    char potenciale[] = "potenciale.txt";
    archivo_potenciale = fopen(potenciale, "w");

    //abrimos fichero energia total
    FILE *archivo_;
    char conservacion[] = "energianplanet.txt";
    archivo_ = fopen(conservacion, "w");

     //abrimos fichero L
 FILE *archivo__;
    char momentoangular[] = "angularnplanet.txt";
    archivo__ = fopen(momentoangular, "w");

    //abrimos archivo excentricidad y periodo
    FILE *archivo___;
    char excperiod[] = "excperiodnplanet.txt";
    archivo___ = fopen(excperiod, "w");

  //definicion de variables
    int i, step, filas, j, k, N;
    long double *m, *r_x, *r_y, *v_x, *v_y,*w_x, *w_y ,*a_x, *a_y, *period, *periodref,*exc,energiatotal, momentoang, g, auxp, errorr, incd, V, cin;
    double h, t;

      //definimos el paso
    h=0.01;
    t=1000;
    step=t/h;

 //numero de planetas filas=N+1 donde N es nº planetas ojo cuenta con el sol
    N=2; //numero de planetas que se le suma a la simulacion
    filas=10 + N;
    //con pluto filas=10;
   
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



// distancia que van a tener entre los N planetas
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
//le ponemos la posicion al resto de N planetas
    for(i=10;i<filas;i++){
    r_x[i]=r_x[i-1] + incd;
    }

//reescalamos r_x
  rescr(r_x, filas);

// m: Masas de los planetas en kilogramos
m[0] = masaSOL;
m[1] = masaMER;
m[2] = masaVEN;
m[3] = masaTIE;
m[4] = masaMAR;
m[5] = masaJUP;
m[6] = masaSATURN;
m[7] = masaURANUS;
m[8] = masaNEPTUNE;
m[9] = masaPLUTO;
//le ponemos masa al resto de N planetas
  for(i=10;i<filas;i++){
    m[i]= masaJUP; //en este caso, la masa de Jupiter
    }

//reescalamos masas
   rescm(m, filas);

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
//le ponemos velocidad inicial al resto de N planetas
  for(i=10;i<filas;i++){
    v_y[i]= VPLUTO; //en este caso, la velocidad de pluton
    }

//reescalamos la velocidad como v=x/t reescalar v es dividir por la posicion y multiplicar por el reesc de t
rescr(v_y, filas);
corriget(v_y, filas);

    //inicializar la componente "y" y el vector "a", "w", "period" y "v" a 0
    for(i=0; i<filas; i++)
    {
        r_y[i]=0;
        v_x[i]=0;
        a_x[i]=0;
        a_y[i]=0; 
        w_x[i]=0;
        w_y[i]=0;
     }
    
    //llamar a la aceleracion para la primera iteracion
    aceleracion(a_x, m, r_x, r_y, filas);
    aceleracion(a_y, m, r_y, r_x, filas);

    /*ALGORITMO DE VERLET*/
    for(int a=0; a<step; a++)
    {
        
       //calculo nueva posicion x e y
        calculopos(r_x, v_x, a_x, h, filas);
        calculopos(r_y, v_y, a_y, h, filas);
  
        //calculo w
        calculow(w_x,v_x,a_x, h, filas);
        calculow(w_y,v_y,a_y, h, filas);

        //calculo la nueva aceleración para x e y
        aceleracion(a_x, m, r_x, r_y, filas);
        aceleracion(a_y, m, r_y, r_x, filas);

        //calculo la nueva velocidad
        calculov(v_x , w_x, a_x, h, filas);
        calculov(v_y , w_y, a_y, h, filas);

        //reescalo r para escribirlo en el fichero
        corriger(r_x, filas);
        corriger(r_y, filas);


         //imprimo posiciones nuevas
         for(k=0; k<filas; k++)
        {
            fprintf(archivo, "%.10Lf", r_x[k]);
            fprintf(archivo, ", %.10Lf\n", r_y[k]);
        }
        fprintf(archivo, "\n");

        //reescalamos r de nuevo
        rescr(r_x, filas);
        rescr(r_y, filas);

         //calculamos energia y momento angular total
        V=potencial(m, r_x, r_y, a_x, a_y,filas);
        cin=cinetica(m, v_x, v_y, filas);
        energiatotal=  cin+V;
        momentoang= mangular(m, r_x, r_y,  v_x, v_y, filas);

        //escribimos etotal en el fichero
           fprintf(archivo_, "%.10Lf\n", energiatotal);
        fprintf(archivo_, "\n");

        //escribimos cinetica en el fichero
           fprintf(archivo_cinetic, "%.10Lf\n", cin);
        fprintf(archivo_cinetic, "\n");

        //escribimos potencial en el fichero
           fprintf(archivo_potenciale, "%.10Lf\n", V);
        fprintf(archivo_potenciale, "\n");

        //escribimos momento angular en el fichero
           fprintf(archivo__, "%.10Lf\n", momentoang);
        fprintf(archivo__, "\n");
    }


    // Liberar la memoria asignada para los vectores
    free(m);
    free(r_x);
    free(r_y);
    free(v_x);
    free(v_y);
    free(w_x);
    free(w_y);
    free(a_x);
    free(a_y);

    //cerramos ficheros
 fclose(archivo);
 fclose(archivo_);
 fclose(archivo__);
 fclose(archivo___);
 fclose(archivo_cinetic);
 fclose(archivo_potenciale);

    return 0;
}

//funcion reescalar masa
void rescm(long double *masa, int n)
 {
    int i;
    for (i = 0; i < n; i++) {
        masa[i] = masa[i] / masaSOL;
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
//funcion que devuelve la escala original de r
void corriger(long double *pos, int n)
{
    for (int i = 0; i < n; i++)
    {
        pos[i] = pos[i] * c;
    }
}

//funcion que reescala tiempo
void reesct(long double *tiempo, int n)
{
    for (int i = 0; i < n; i++)
    {
        tiempo[i] = tiempo[i] * sqrt((G * masaSOL) / pow(c, 3));
    }
}

//funcion que devuelve la escala original de tiempo
void corriget(long double *tiempo, int n)
{
    for (int i = 0; i < n; i++)
    {
        tiempo[i] = tiempo[i] / sqrt((G * masaSOL) / pow(c, 3));
    }
}

//funcion aceleracion imput; masa, posicion de la aceleracion que se va a calcular, la otra posicon, filas. output; aceleracion que se va a calcular
void aceleracion(long double *aceleracion, long double *masa, long double *posx, long double *posy, int n)
{

    //F(r)= -parcialV/r = 24eps*sigma^6(2sigma^6-r^6)/r^13
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

//funcion que calcula el momento angular total
long double mangular(long double *masa, long double *posx, long double *posy, long double *velx, long double *vely, int n)
{
long double L;
L=0;
     for (int i = 1; i < n; i++)
    {
        L += masa[i] * sqrt(pow(posx[i], 2) + pow(posy[i], 2)) * sqrt(pow(velx[i], 2) + pow(vely[i], 2));
    }
    return L;
}

//funcion que calcula energia cinetica total
long double cinetica(long double *masa, long double *velx, long double *vely, int n)
{
    long double kinetic=0;

    for(int i=1; i<n; i++)
    {
     kinetic += 0.5 * masa[i] * (velx[i]*velx[i] + vely[i]*vely[i]);
    }

    return kinetic;
}

//funcion que calcula la energia potencial total
long double potencial(long double *masa, long double *posx, long double *posy, long double *acelx, long double *acely,int n)
{
   long double pot=0;

   //iniciamosvector cinetica y potencial
    for(int i=1; i<n; i++)
    {
     pot+= -masa[i] * sqrt(pow(posx[i], 2) + pow(posy[i], 2)) * sqrt(pow(acelx[i], 2) + pow(acely[i], 2));
    }

    return pot;
}
