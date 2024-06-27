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
    int i, N, j, n, nciclos, p, m, mT;
    double complex *fonda, *beta, *alpha, *A_0, *chi;
    double lambda, h, x_0, sigma, *V, k_0, s, norm, numaleatorio, evaluado, posiblet, anteriort;

     //Inicializo el valor de la serie de números aleatorios
    srand(time(NULL));

    //abrimos arhcivo de salida
    FILE *archivo;
    archivo = fopen("sch.dat", "w");

    if (archivo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }\
    //abrimos arhcivo norma
    FILE *normadata;
    normadata = fopen("norma.dat", "w");

    if (normadata == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    
    //definimos N
    N = 1000/2;
    p=1000;

    // Asignar memoria para los vectores
    alpha = (double complex *)malloc((N + 1)*sizeof(double complex));
    V = (double *)malloc((N + 1)*sizeof(double));
    fonda = (double complex *)malloc((N + 1)*sizeof(double complex));
    chi = (double complex *)malloc((N + 1)*sizeof(double complex));
    beta =(double complex *)malloc((N + 1)*sizeof(double complex));


    //parametros a cambiar: lambda, nciclos y N

    // Definicion de constantes
    //nciclos=50.;
    nciclos=50.;
    h=0.01;
    norm=0.;
    lambda = 0.5;
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
m=100;
for(int expe=0; expe<=m; expe++)
{
    //inicalizamos los valores que van a comparar
    anteriort=0.;
    posiblet=0.1;

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

        // 5. n = n + 1, ir a al paso 2
         //imprimimos la funcion de onda
            for(j=0;j<=N;j++){
             //   fprintf(archivo, "%i, %f, %f, %f, %f\n",j, cabs(fonda[j]), creal(fonda[j]),cimag(fonda[j]), V[j]); //todos paramretros
                fprintf(archivo, "%i, %f, %f\n",j, cabs(fonda[j]), V[j]); //solo vabs
            }
             fprintf(archivo, "\n");

             // Calculamos la probabilidad a la derecha (primer máximo local t=742) (cuando el anterior deje de ser menor que el nuevo se llega al máximo local)
            posiblet = detectorD(fonda, N);

            if(posiblet < anteriort)
            {
                n = p+1;
            }
            else
            {
                anteriort = posiblet;
            }

        }

       // 6.1 Simulamos el proceso de medicion generando un numero aleatorio p ∈ [0, 1]. Si p > PD(nD) habremos detectado 
       //la particula y actualizamos el valor mT = mT + 1. Si p < PD(nD) no se habria detectado la particula y
        //actualizamos mT = mT + 0.
        numaleatorio=aleatorio();

        if (numaleatorio>posiblet) {
            mT++;
        }
        else mT=mT+0;

}
printf("%f", mT/m); 

    //cerramos fichero de salida
        fclose(archivo);
        fclose(normadata);

    // Liberar la memoria dinamica
    free(V);
    free(alpha);
    free(chi);
    free(fonda);
    free(beta);

    return 0;
}

double detectorD(double complex *psi, int N)
{
    int j;
    double PD=0.;

    for (j=4.*N/5; j<=N; j++)
    {
        PD+=pow(cabs(psi[j]),2);
    }

    return  PD;
}

double aleatorio() 
{
    return (double)rand() / (double)RAND_MAX;
}