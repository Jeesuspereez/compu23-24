#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define PI (3.1415926535897932384)

double detectorD(double complex *psi, int N);
double aleatorio();

int main(void)
{
    //declaracion de variables
    int i, N, j, n, nciclos, p, m, mT, chivato, chivato2, y;
    double complex *fonda, *beta, *alpha, *A_0, *chi;
    double lambda, h, x_0, sigma, *V, k_0, s, norm, numaleatorio, evaluado, posiblet, anteriort, aux, *tiempillo, K, apartadotres;

     //Inicializo el valor de la serie de números aleatorios
    srand(time(NULL));

    //chivato
    chivato=0;
    chivato2=0;
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

    //abrimos arhcivo de salida
    FILE *prueba;
    prueba = fopen("tiempillo.txt", "w");

    if (prueba == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    
    //definimos N
    N = 2000;
  //  p=1000;
  p=10000;

    // Asignar memoria para los vectores
    alpha = (double complex *)malloc((N + 1)*sizeof(double complex));
    V = (double *)malloc((N + 1)*sizeof(double));
    tiempillo = (double *)malloc((p + 10)*sizeof(double));
    fonda = (double complex *)malloc((N + 1)*sizeof(double complex));
    chi = (double complex *)malloc((N + 1)*sizeof(double complex));
    beta =(double complex *)malloc((N + 1)*sizeof(double complex));

    //inicializamos tiempillo
    for( y=0; y<p+1 ; y++)
    {
        tiempillo[y]=0.;
    }


    //parametros a cambiar: lambda, nciclos y N

    // Definicion de constantes
    //nciclos=50.;
    nciclos=50.;
    h=0.01;
    norm=0.;
    lambda = 0.05;
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

    //hacemos el experimento m veces
m=1000;
apartadotres=0.;
for(int expe=0; expe<m; expe++)
{

        for(y=0; y<p+1 ; y++)
    {
        tiempillo[y]=0.;
    }
    chivato=0;

        for(n=0; n<p+1; n++){ 
       // tiempillo[n+1]=0.;
    //    tiempillo[0]=0.; // lo acabo de quitar pero creo q no hace falta

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

            // 5. n = n + 1, ir a al paso 2
            //imprimimos la funcion de onda
             for(j=0;j<=N;j++){
             //   fprintf(archivo, "%i, %f, %f, %f, %f\n",j, cabs(fonda[j]), creal(fonda[j]),cimag(fonda[j]), V[j]); //todos paramretros
                fprintf(archivo, "%i, %f, %f\n",j, cabs(fonda[j]), V[j]); //solo vabs
            }
             fprintf(archivo, "\n");

             // Calculamos la probabilidad a la derecha
            tiempillo[n+1]=detectorD(fonda,N);

            if(tiempillo[n+1]>=tiempillo[n])
            {
                aux=tiempillo[n+1];
                chivato++;
            }

            else n=p+1;

            tiempillo[n]=aux;
            chivato2++;

        /*   if(expe==10){
                  for(j=0;j<=p+9;j++)
                  {
                fprintf(prueba, "%f,", tiempillo[j]);
            }

            }
            */ 

        }

       // 6.1 Simulamos el proceso de medicion generando un numero aleatorio p ∈ [0, 1]. Si p > PD(nD) habremos detectado 
       //la particula y actualizamos el valor mT = mT + 1. Si p < PD(nD) no se habria detectado la particula y
        //actualizamos mT = mT + 0.
      numaleatorio=aleatorio();

    /*   if (numaleatorio>tiempillo[chivato]) {
            mT++;
        }

        */ 
       apartadotres+=tiempillo[chivato];

        fprintf(prueba, "%i, %f", chivato, tiempillo[chivato]);
        fprintf(prueba, "\n");

       if (numaleatorio<tiempillo[chivato]) {
            mT++;
        }

 //       printf("%.20f,", tiempillo[chivato]); 

}

K=(0.0+mT)/m;
apartadotres=apartadotres/m;
printf("%f,", K); 
// printf("%i,", chivato);
printf("%i,", chivato2);
printf("%i,", mT);
printf("%i,", m);
printf("\n"); 
printf("apartado tres (promedio de Pd(nd)): ");
printf("%f", apartadotres); 

    //cerramos fichero de salida
        fclose(archivo);
        fclose(normadata);
        fclose(prueba);

    // Liberar la memoria dinamica
    free(V);
    free(alpha);
    free(chi);
    free(fonda);
    free(beta);
    free(tiempillo);

    return 0;
}

double detectorD(double complex *psi, int N)
{
    int j;
    double PD=0.;

    for (j=4*N/5; j<=N; j++)
    {
        PD+= cabs(psi[j])*cabs(psi[j]);
    }

    return  PD;
}

double aleatorio() 
{
    return (double)rand() / (double)RAND_MAX;
}