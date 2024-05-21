#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define PI (3.1415926535897932384)

int main(void)
{
    //declaracion de variables
    int i, N, j, nciclos, n, ciclos;
    double complex *fonda, *beta, *alpha, *A_0, *chi;
    double lambda, h, x_0, sigma, *V, k_0, s;

    ciclos=1000;

    //definimos N
    N = 1000;

    //definimos nciclos
    nciclos = N/16;

     //definimos h
    h=0.01;

    FILE *archivo;
    archivo = fopen("sch.dat", "w");

    if (archivo == NULL) {
        printf("Error al abrir el archivo.\n");
        return 1;
    }

    // Asignar memoria para los vectores
    fonda = (double complex *)malloc((N+1) * sizeof(double complex));
    beta = (double complex *)malloc((N+1) * sizeof(double complex));
    V = (double *)malloc((N+1) * sizeof(double));
    alpha = (double complex *)malloc((N+1) * sizeof(double complex));
    A_0 = (double complex *)malloc((N+1) * sizeof(double complex));
    chi = (double complex *)malloc((N+1) * sizeof(double complex));


    // Definimos constantes
    lambda = 0.1;
    k_0 = (2.*PI*nciclos)/N;
    s = 1./(4.*k_0 *k_0);
    x_0 = N*h/4.;
    sigma = N*h/16.;

    // Condición de contorno para fonda
    fonda[0] =0.;
    fonda[N] = 0.;

    alpha[N-1] = 0.0;
    beta[N-1] = 0.0;
    chi[N] = 0.0;
    chi[0] = 0.0;

    // Condiciones de contorno para fonda
    for(j = 1; j < N-1; j++)
    {
        fonda[j] = cexp(k_0 * I * j +0.) * cexp(-1.*(j - x_0) * (j - x_0) / (2. * sigma * sigma));
    }

    // Potencial Vj vale 0 fuera del intervalo 2N/5,3N/5
    for(i = 0; i < 2*(N/5); i++)
    {
        V[i] = 0.;
    }

    for(i = 2*(N/5); i < 3*N/5; i++)
    {
        V[i] = lambda * k_0 * k_0;
    }

    for(i = 3*(N/5); i < N; i++)
    {
        V[i] = 0.;
    }
    

   //calcular V_j
   /* for(j=0;j<N;j++){
        if((j>=2*N/5)&&(j<=3*N/5)){
            V[j]=1.0*lambda*k_0*k_0;
        }
        else{
            V[j]=0.0;
        }

    }
*/
    alpha[N-2]=0.;
    // Calculamos una vez sola alpha ya que no depende del tiempo
    for ( j = N-2; j > 0; j--)
    { 
        alpha[j-1] =-1.0 / ((-2.0 + (2.0 * I / s) - V[j]) + alpha[j]);
    }

        for(n=0; n<1000; n++){

            //imprimimos la funcion de onda
            for(j=0;j<N;j++){
                fprintf(archivo, "%d, %f, %f, %f, %f\n",j, creal(fonda[j]),cimag(fonda[j]), cabs(fonda[j]), V[j]);
         /*   fprintf(archivo, "%d", j);
            fprintf(archivo, ", %.10f", creal(fonda[j]));
            fprintf(archivo, ", %.10f", cimag(fonda[j]));
            fprintf(archivo, ", %.10f", cabs(fonda[j]));
            fprintf(archivo, ", %.10f\n", V[j]);
            */
            }

        //contorno para alpha y beta
       // alpha[N-1]=beta[N-1]=0.;  

        // 2. Calcular beta utilizando la recurrencia (22).
          //  beta[N-2]=0.;
            for(j=N-2;j>0; j--){
            A_0[j] = -2. - V[i] + 2. * I / s;
            beta[j-1]=(4.0*i*(fonda[j]/s)-beta[j])/(A_0[j]+alpha[j]);
            }

            // 3. Calcular chi a partir de (20).
                chi[0]=0.;    
                for(j=0;j<N-2;j++){
                chi[j+1]= alpha[j]*chi[j]+beta[j];
                 }  

        // 4. Calcular phi(j,n+1) de (15).
        for(j=0;j<N;j++){
            fonda[j]=chi[j]-fonda[j];
        }

        // 5. n = n + 1, ir a al paso 2

        }

        fclose(archivo);

    // Liberar la memoria
    free(fonda);
    free(beta);
    free(V);
    free(alpha);
    free(A_0);
    free(chi);

    return 0;
}
