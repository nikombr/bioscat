#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "../../lib/GaussianProcess.h"
#include <cuda_runtime_api.h>
#include <ctime>
extern "C" {

void print_matrix(double **M_h, double * M_log, bool device, int n) {
    if (n < 21) {
        if (device) {
            cudaMemcpy(*M_h, M_log, n * n * sizeof(double), cudaMemcpyDeviceToHost);
        }
        printf("\n");
        printf("<");
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                if (i != n-1) printf("%.8f, ",M_h[k][i]);
                else printf("%.8f",M_h[k][i]);
            }
            if (k != n-1) printf(";\n");
            
        }
        printf(">");
        printf("\n");
        printf("\n");
    }

}


void profileGaussianProcess(double * x, double * y, int n, double * hyper, int num, int dim, int dev, int type_covfunc) {
    double start, stop;

    FILE *file;
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");
    char filename[256];
    if (type_covfunc == 1) sprintf(filename, "%s/Results/profiling/GaussianProcess_squared_exponential.txt",dir);
    else sprintf(filename, "%s/Results/profiling/GaussianProcess_matern.txt",dir);
    file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    int iter = 1e4;

    fprintf(file,"%d ", n);

    start = omp_get_wtime();
    GaussianProcess GP = GaussianProcess(x, y, n, hyper, num, dim, dev, type_covfunc);
    stop = omp_get_wtime();

    fprintf(file,"%.8e ", stop - start);

    start = omp_get_wtime();
    GP.covariance_matrix();
    stop = omp_get_wtime();

    fprintf(file,"%.8e ", stop - start);

    start = omp_get_wtime();
    GP.cholesky();
    stop = omp_get_wtime();

    fprintf(file,"%.8e ", stop - start);

    start = omp_get_wtime();
    // Seed the random number generator with the current time
    srand(time(NULL));
    //srand(0);
    for (int k = 0; k < iter; k++) {
        GP.realisation();
    }
    stop = omp_get_wtime();

    fprintf(file,"%.8e ", (stop - start)/iter);

    GP.free();

    fprintf(file,"\n");


}

}

