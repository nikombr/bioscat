#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
extern "C" {
#include "../lib/GaussianProcess.h"
#include <cuda_runtime_api.h>

void gaussian_process_inner(double * x, double * y, int n, double * hyper, int num, int dim) {
    double start, stop;


    start = omp_get_wtime();
    GaussianProcess GP = GaussianProcess(x, y, n, hyper, num, dim);
    stop = omp_get_wtime();

    printf("Initialization and allocation: %.4f seconds\n", stop - start);

    start = omp_get_wtime();
    GP.covariance_matrix();
    stop = omp_get_wtime();

    printf("Computing covariance matrix: %.4f seconds\n", stop - start);
    if (n < 21) {
        if (GP.device) {
            cudaMemcpy(*GP.M_h, GP.M_log, n * n * sizeof(double), cudaMemcpyDeviceToHost);
        }
        printf("\n");
        printf("<");
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                if (i != n-1) printf("%.8f, ",GP.M_h[k][i]);
                else printf("%.8f",GP.M_h[k][i]);
            }
            if (k != n-1) printf(";\n");
            
        }
        printf(">");
        printf("\n");
        printf("\n");
    }

    start = omp_get_wtime();
    GP.cholesky();
    stop = omp_get_wtime();

    printf("Cholesky factorization: %.4f seconds\n", stop - start);

}

}