#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <omp.h>

#define PI 3.14159265359
#define G (long double)(6.674 * pow(10, -11)) // Constante de gravitación universal en m^3 kg^-1 s^-2
#define masaTIE (long double)(5.972 * pow(10, 24)) // Masa de la tierra en kg
#define masaLUN (long double)(7.349 * pow(10, 22)) // Masa de la luna en kg
#define masaCOHETE (long double)(1000) // Masa del cohete en kg
#define OMEGA (long double )(2.6617*pow(10,-6)) //velocidad angular de la luna en rad/s
#define radioTIERRA (long double )(6.37816*pow(10,6)) //radio de la tierra en kg
#define radioLUNA (long double )(1.7374*pow(10,6)) // Radio de la luna en kg
#define distanciaTIELUN (long double )(3.844*pow(10,8)) // Distancia entre la tierra y la luna 

long double f1(long double P_r);
long double f2(long double r, long double p_phi);
long double f3(long double r, long double p_phi, long double phi, long double w, long double nu, long double delta, long double t);
long double f4(long double r, long double phi, long double w, long double nu, long double delta, long double t);

int main(void)
{
    int i, filas,j;
    filas=5;
    long double **k,phi, r, h, Pr, Pphi, mu, delta,t, Vesc,v, velx,vely, maximo, rdot, conservacionprog, angulo;

    j=0;
    t=0.;
    conservacionprog=0.;

     FILE *archivo;
    char rock[] = "cohete.txt";
    archivo = fopen(rock, "w");

  /*  FILE *archivo_;
    char conservacion[] = "conservacion.txt";
    archivo_ = fopen(conservacion, "w"); */

    //definimos la matriz de k, la fila indica k1, k2 k3 k4 mientras que la columna indica la funcion
    k = (long double **)malloc((filas+1) * sizeof(long double *));
    if (k == NULL) {
        printf("Error: no se pudo asignar memoria para las filas\n");
        return 1;
    }

    // Asignar memoria para las columnas de cada fila
    for (i = 0; i < filas+1; i++) {
        k[i] = (long double *)malloc((filas+1) * sizeof(long double));
        if (k[i] == NULL) {
            printf("Error: no se pudo asignar memoria para las columnas de la fila %d\n", i);
            return 1;
        }
    }

    //dar valor a constantes
    h=0.01;
	delta = (1.*G * masaTIE) / pow(distanciaTIELUN, 3.0);
	mu = 1.*masaLUN / masaTIE;
    //maximo=10000000;
    maximo=9*10000000;


	//Condiciones del cohete y reescalamiento
	r = 1.0*radioTIERRA / distanciaTIELUN;
	//v=1.0*11160.0/distanciaTIELUN;
    v=1.0*11115.0/distanciaTIELUN; //cambiar la velocidad aqui

	angulo = PI/6.0;  //cambiar el angulo
    phi=0.;
	velx = v * cosl(PI / 4.0); // PI/4.0 es la posición desde la superficie a la Luna
	vely = v * sinl(PI / 4.0); // PI/4.0 es la posición desde la superficie a la Luna
	Pr = v * cosl(angulo);
	Pphi = r * v * sinl(angulo);
	rdot = sqrt(1 + r * r - 2.0 * r * cosl(phi- OMEGA*0.));

     fprintf(archivo, "%Lf", 1.0*r*cosl(phi));
              fprintf(archivo, ", %Lf\n", 1.0*r*sinl(phi));
               fprintf(archivo, "%Lf", 1.0*cosl(OMEGA*t*1.0));
                fprintf(archivo, ", %Lf\n", 1.0*sinl(OMEGA*t*1.0));
                fprintf(archivo, "\n");

   for(j=0;j<maximo;j++)
    {
        if(j%5000==0){
             fprintf(archivo, "%Lf", 1.0*r*cosl(phi));
              fprintf(archivo, ", %Lf\n", 1.0*r*sinl(phi));
               fprintf(archivo, "%Lf", 1.0*cosl(OMEGA*t*1.0));
                fprintf(archivo, ", %Lf\n", 1.0*sinl(OMEGA*t*1.0));
                 fprintf(archivo, "\n");
           //      fprintf(archivo_, " %Lf,",conservacionprog); PARA VER CONSERVACION QUITAR RECUERDA QUITAR EL COMENTARIO DE FCLOSE Y LA DECLARACION   
             //    fprintf(archivo_, "\n");
        }

    //2. Evaluar k(1) 
    k[1][1]=h*f1(Pr);
    k[1][2]=h*f2(r, Pphi);
    k[1][3]=h*f3(r, Pphi, phi, OMEGA, mu, delta, t);
    k[1][4]=h*f4(r, phi, OMEGA, mu, delta, t); 

    //3. Evaluar k(2)
    k[2][1]=h*f1(Pr+ k[1][3] / 2.);
    k[2][2]=h*f2(r + k[1][1] / 2.0,Pphi + k[1][4] / 2.0);
    k[2][3]=h*f3(r + k[1][3] / 2.0,Pphi + k[1][4] / 2.0, phi + k[1][2] / 2.0, OMEGA, mu, delta, 1.0*j*h);
    k[2][4]=h*f4( r + k[1][1] / 2.0, phi + k[1][2] / 2.0, OMEGA, mu,delta, 1.0*j*h);

    //4. Evaluar k(3)
    k[3][1]=h*f1(Pr+ k[2][3] / 2.);
    k[3][2]=h*f2(r + k[2][1] / 2.0,Pphi + k[2][4] / 2.0);
    k[3][3]=h*f3(r + k[2][3] / 2.0,Pphi + k[2][4] / 2.0, phi + k[2][2] / 2.0, OMEGA, mu, delta, 1.0*j*h);
    k[3][4]=h*f4( r + k[2][1] / 2.0, phi + k[2][2] / 2.0, OMEGA, mu,delta, 1.0*j*h);

    //5. Evaluar k(4)
    k[4][1]=h*f1(Pr+ k[3][3]);
    k[4][2]=h*f2(r + k[3][1],Pphi + k[3][4]);
    k[4][3]=h*f3(r + k[2][1],Pphi + k[3][4] , phi + k[3][2], OMEGA, mu, delta, j*h);
    k[4][4]=h*f4( r + k[3][1], phi + k[3][2], OMEGA, mu,delta, j*h);

    //6. yn(t+h)
        Pphi= Pphi + (k[1][4]+2.0*k[2][4]+2.0*k[3][4]+k[4][4])/6.0;
        r=r + (k[1][1]+2.0*k[2][1]+2.0*k[3][1]+k[4][1])/6.0;
        phi= phi + (k[1][2]+2.0*k[2][2]+2.0*k[3][2]+k[4][2])/6.0;
        Pr= Pr + (k[1][3]+2.0*k[2][3]+2.0*k[3][3]+k[4][3])/6.0;

        //Vemos si se conserva
        conservacionprog=(Pr*Pr/2.0) + (Pphi*Pphi/(2.0*r*r)) - ((G/(pow(distanciaTIELUN,3)))*((masaTIE/r) + (masaLUN/(sqrt(1.0 + r*r - 2.0*r*cos(phi - 1.0*OMEGA*t))))));

    //7. t=t+h
        t=t+h;
    }
    



// Liberar la memoria
    for (int i = 0; i < filas+1; i++) {
        free(k[i]);
    }
    free(k);

    return 0;
        fclose(archivo);
        // fclose(archivo_);

}

//funcion que devuelve r derivado
long double f1(long double P_r)
{
    return P_r;
}

//funcion que devuelve phi derivado
long double f2(long double r, long double p_phi)
{
    long double phi;

    phi=1.0*p_phi/(r*r);

    return phi;
}

//funcion que devuelve P_r derivado
long double f3(long double r, long double p_phi, long double phi, long double w, long double nu, long double delta, long double t)
{
    long double P_r;
    long double aux;

    aux=sqrt(1+r*r -2.*r*cos(phi - 1.0*w*t));
    P_r= (1.0*p_phi*p_phi)/(pow(r,3)) -delta*(1./(r*r) + nu/pow(aux,3)*(r -cos(phi -1.0*w*t)) );

    return P_r;
}


//funcion que devuelve P_phi derivado
long double f4(long double r, long double phi, long double w, long double nu, long double delta, long double t)
{
    long double P_phi;
    long double aux;

    aux=sqrt(1.+r*r -2.*r*cos(phi - 1.0*w*t));
    P_phi= (-1.0*delta*nu*r*sin(phi-1.0*w*t))/pow(aux,3);

    return P_phi;
}

