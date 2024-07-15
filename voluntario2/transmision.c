#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>
#include <stdbool.h>

#define PI (3.1415926535897932384)

double detectorD(double complex *psi, int N);
double aleatorio();

int main(void)
{
    //declaracion de variables
    int i, N, j, n, nciclos, p, m, mT, y, zz;
    double complex *fonda, *beta, *alpha, *A_0, *chi, *operadorenergia, *fondaaux, *derivadafonda, *segundaderiv;
    double lambda, h, x_0, sigma, *V, k_0, s, norm, numaleatorio, aux, K, apartadotres, pd, auxiliar, normamedir, posicion, cinetica, iteracion, normapos;
    int find, activarpos;
    double iteracioncinetica=0;
    double normacinetica=0;

    //inicializamos valores
    find=0;
    pd=0;
    posicion=0.;
    cinetica=0;

    //si queremos calcular la posicion ponemos activarpos=1 si no le ponemos cualq otro valor por ejemplo 0
    activarpos=1;

     //Inicializo el valor de la serie de números aleatorios
    srand(time(NULL));

    //abrimos arhcivo de salida
    FILE *archivo;
    archivo = fopen("sch.dat", "w");

    if (archivo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }
    //abrimos arhcivo norma
    FILE *normadata;
    normadata = fopen("norma.dat", "w");

    if (normadata == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    FILE *pos;
    pos = fopen("posicion.txt", "w");

    if (pos == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    //abrimos arhcivo de salida
    FILE *energy;
    energy = fopen("cinetica.txt", "w");

    if (energy == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }


    //abrimos arhcivo de salida
    FILE *tiempo;
    tiempo = fopen("tiempo.txt", "w");

    if (tiempo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    
    //definimos N
    N = 1000;
    //  p=1000;
    p=2000;

    // Asignar memoria para los vectores
    alpha = (double complex *)malloc((N + 1)*sizeof(double complex));
    V = (double *)malloc((N + 1)*sizeof(double));
    fonda = (double complex *)malloc((N + 1)*sizeof(double complex));
    chi = (double complex *)malloc((N + 1)*sizeof(double complex));
    beta =(double complex *)malloc((N + 1)*sizeof(double complex));
    operadorenergia = (double complex *)malloc((N + 1)*sizeof(double complex));
    fondaaux = (double complex *)malloc((N + 1)*sizeof(double complex));
    derivadafonda = (double complex *)malloc((N + 1)*sizeof(double complex));
    segundaderiv = (double complex *)malloc((N + 1)*sizeof(double complex));

    //parametros a cambiar: lambda, nciclos y N

    // Definicion de constantes
    //nciclos=50.;
    nciclos=50.;
    h=0.01;
    norm=0.;
    lambda = 0.4;
    k_0 = (2.*PI*nciclos)/(N+0.);
    s = 1./(4.*k_0 *k_0);
    //parametros distribucion
    x_0 = N/8.;
    sigma = N/20.;

    // Condición de contorno para la funcion de onda
    fonda[0] =0.;
    fonda[N] = 0.;

    // Condiciones de contorno para fonda
    for(j = 1; j < N; j++)
    {
        fonda[j] = cexp(1.0 * I * k_0 * j) * cexp(-1.0 * (j - x_0) * (j - x_0) / (2.0 * sigma * sigma));
         norm= norm+ cabs(fonda[j])*cabs(fonda[j]);
    }

    //normalizamos la funcion de onda
     for(j = 1; j < N; j++){
        fonda[j] = fonda[j]/sqrt(norm);
    }

    // Potencial Vj vale 0 fuera del intervalo 2N/5,3N/5
    for(i = 0; i < (2 * N / 5); i++)
    {
        V[i] = 0.;
    }

    for(i = (2 * N / 5); i <= (3 * N / 5) ; i++)
    {
        V[i] = k_0 * k_0* lambda;
    }

    for(i = (3 * N / 5) + 1; i <= N; i++)
    {
        V[i] = 0.;
    }

    //condiciones de contorno para los auxiliares
    chi[N] = 0.0;
    chi[0] = 0.0;
    alpha[N-1] = 0.0;
    beta[N-1] = 0.0;

    // Calculamos una vez sola alpha ya que no depende del tiempo
    for ( j = N-2; j > 0; j--)
    { 
        alpha[j-1] = -1.0 / ((-2.0 + (2.0 * I / s) - V[j]) + alpha[j]);
    }

        //calculamos la funcion de onda
        for(n=0; n<p+1; n++){

        // 2. Calcular beta utilizando la recurrencia (22).
            for(j=N-2;j>0; j--){
           beta[j-1] = (1.0 / (-2.0 + 2.0 * (I/s) + alpha[j] - V[j]))*(4.0 * I * (fonda[j]/s) - beta[j]);
            }

        // 3. Calcular chi a partir de (20).  
            for(j=0;j<=N-2;j++){
                chi[j+1]= beta[j] + alpha[j]*chi[j];
            }  

        // 4. Calcular phi(j,n+1) de (15).
        norm=0.;
            for(j=0;j<=N; j++){
                fonda[j]=chi[j]-fonda[j];
                norm += cabs(fonda[j])*cabs(fonda[j]);
            }
            fprintf(normadata, "%f", norm);
            fprintf(normadata, "\n");

            //calculamos posicion y energia

            //posicion
                iteracion=normapos=0.;

                //calculamos el operador
                 for (int g = 0; g <= N; g++)
                {
                    iteracion += cabs(fonda[g])*g;
                    normapos += cabs(fonda[g]);
                }

                posicion=1.0*iteracion/normapos;

            fprintf(pos, "%f\n", posicion);

            //archivo de tiempo para graficar pos vs tiempo y energia vs tiempo
            fprintf(tiempo, "%i\n", n);

            //energia
                //calculamos el operador energia
                for (zz = 0; zz <= N; zz++)
                    {
                        operadorenergia[zz] =creal(fonda[zz]) - I*cimag(fonda[zz]);
                    }

                //calculamos la derivada 
                for (int i = 0; i <= N; i++)
                {
                    if (i!=N)
                    {
                        derivadafonda[i] = fonda[i+1] - fonda[i];
                    }
                    else
                    {
                        derivadafonda[i] = 0.0;
                    }
        
                }

                //y su segunda derivada
                for (int g = 0; g <= N; g++)
                {
                    if (g!=N)
                    {
                        segundaderiv[g] = derivadafonda[g+1] - derivadafonda[g]; 
                    }
                    else
                    {
                        segundaderiv[g] = 0.0;
                    }
        
                }

            iteracioncinetica=0;
            normacinetica=0;
            //calculamos definitivamente la energia cinetica
            for (int i = 0; i <= N; i++)
            {
                iteracioncinetica += -creal(operadorenergia[i]*segundaderiv[i]);
                normacinetica += cabs(operadorenergia[i]);
            }

            cinetica = iteracioncinetica/(2*normacinetica); 

            fprintf(energy, "%f\n", cinetica);

            //buscamos el maximo
            //calculamos la norma total
            normamedir = 0.0;

            for (int i = 0; i <= N ; i++) 
            {
            normamedir += cabs(fonda[i]) * cabs(fonda[i]);
            }

            auxiliar = detectorD(fonda, N);
        
            if (find == 0 &&  auxiliar >= pd) 
            {
                pd = auxiliar;

             //   t = j;
            } 
            else if ((pd - auxiliar) > 0.001) 
            {
                find = 1;
            } 

        }

            pd = pd / normamedir;
    
            mT = 0;

            m=1000;
            //proceso de medicion
             for (int k = 0; k < m; k++) 
            {
                numaleatorio = aleatorio(); //Genero un número aleatorio
        
                if (numaleatorio < pd) 
                {
                    mT = mT + 1;
                }
            }

K=(0.0+mT)/m;
printf("para lambda ");
printf("%f,", lambda); 
printf(" el coeficiente calculado es: ");
printf("%.10f,", K); 
printf("\n"); 

printf("%i,", mT);

printf("apartado tres (promedio de Pd(nd)): ");
printf("%f", pd); 

    //cerramos fichero de salida
        fclose(archivo);
        fclose(normadata);
        fclose(energy);
        fclose(pos);
        fclose(tiempo);

    // Liberar la memoria dinamica
    free(V);
    free(alpha);
    free(chi);
    free(fonda);
    free(beta);
    free(operadorenergia);
    free(derivadafonda);
    free(fondaaux);
    free(segundaderiv);

    return 0;
}

double aleatorio() 
{
    return (double)rand() / (double)RAND_MAX;
}

double detectorD(double complex* psi, int N) 
{
    double PD = 0;
    int i;

    for ( i = N - 100; i <= N ; i++) 
    {
        PD += cabs(psi[i]) * cabs(psi[i]);
    }

    return PD;
}