#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
extern "C" {
#include "../lib/GaussianProcess.h"



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

    start = omp_get_wtime();
    GP.cholesky();
    stop = omp_get_wtime();

    printf("Cholesky factorization: %.4f seconds\n", stop - start);

}

}