#include <stdio.h>
#include <math.h>
#include <stdlib.h>

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
    int i, filas;
    filas=4;
    long double **k,phi, r, h, Pr, Pphi, mu, delta,t, Vesc,v, velx,vely, maximo;

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
	delta = (G * masaTIE) / pow(distanciaTIELUN, 3.0);
	mu = masaLUN / masaTIE;
    maximo=10000;

	//Condiciones del cohete y reescalamiento
	r = radioTIERRA / distanciaTIELUN;
	v = sqrt((2.0*G*masaTIE)/radioTIERRA) / distanciaTIELUN; //reescalada
	phi = 0.49577543479; 
	velx = v * cos(PI / 4.0); // PI/4.0 es la posición desde la superficie a la Luna
	vely = v * sin(PI / 4.0); // PI/4.0 es la posición desde la superficie a la Luna
	Pr = v * cos(PI / 4.0 - phi);
	Pphi = r * v * sin(PI / 4.0 - phi);
	r_ = sqrt(1 + r * r - 2.0 * r * cos(phi)); //fijarse

    while (t<maximo)
    {
       //1. Dar yn = yn(0)

    //2. Evaluar k(1) 
    k[1][1]=h*f1(Pr);
    k[1][2]=h*f2(r, Pphi);
    k[1][3]=h*f3(r, Pphi, phi, OMEGA, mu, delta, t);
    k[1][4]=h*f4(r, phi, OMEGA, mu, delta, t); 

    //3. Evaluar k(2)
    k[2][1]=h*f1(Pr+ k[1][2] / 2.);
    k[2][2]=1;
    k[2][3]=1;
    k[2][4]=1;

    //4. Evaluar k(3)
    k[3][1]=1;
    k[3][2]=1;
    k[3][3]=1;
    k[3][4]=1;

    //5. Evaluar k(4)
    k[4][1]=1;
    k[4][2]=1;
    k[4][3]=1;
    k[4][4]=1;

    //6. yn(t+h)

    //7. t=t+h
        t=t+h;
    }
    
    


// Liberar la memoria
    for (int i = 0; i < filas+1; i++) {
        free(k[i]);
    }
    free(k);

    return 0;

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

    phi=p_phi/(r*r);

    return phi;
}

//funcion que devuelve P_r derivado
long double f3(long double r, long double p_phi, long double phi, long double w, long double nu, long double delta, long double t)
{
    long double P_r;
    long double aux;

    aux=sqrt(1+r*r -2.*r*cos(phi - w*t));
    P_r= (1.0*p_phi*p_phi)/(pow(r,3)) -delta*(1./(r*r) + nu/pow(aux,3)*(r -cos(phi -w*t)) );

    return P_r;
}


//funcion que devuelve P_phi derivado
long double f4(long double r, long double phi, long double w, long double nu, long double delta, long double t)
{
    long double P_phi;
    long double aux;

    aux=sqrt(1.+r*r -2.*r*cos(phi - w*t));
    P_phi= (-1.0*delta*nu*r*sin(phi-w*t))/pow(aux,3);

    return P_phi;
}

