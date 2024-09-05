#include <stdlib.h>
#include <stdio.h>
#include "../lib/cudaMalloc2d.h"
#include <cuda_runtime_api.h>

void gaussian_process(double * x_h, double * y_h, int n, double tau, double ell) {

    // Initialize
    double ** Sigma_h, **L_h, **Sigma_d, **L_d;
    double *Sigma_log, *L_log, *z_h, *z_d, *x_d, *y_d;

    // Allocate matrices
    host_malloc_2d(&Sigma_h, n);
    host_malloc_2d(&L_h, n);
    device_malloc_2d(&Sigma_d, &Sigma_log, n);
    device_malloc_2d(&L_d, &L_log, n);

    // Check allocation
    if (Sigma_h == NULL || L_h == NULL || Sigma_d == NULL || L_d == NULL || Sigma_log == NULL || L_log == NULL) {
        printf("Allocation of matrices failed!");
        return;
    }

    // Allocate vectors
    cudaMalloc((void **) &x_d,     n*sizeof(double));
    cudaMalloc((void **) &y_d,     n*sizeof(double));
    cudaMalloc((void **) &x_d,     n*sizeof(double));
    cudaMallocHost((void **) &z_h, n*sizeof(double));

    // Check allocation
    if (x_d == NULL || y_d == NULL || x_d == NULL || z_h == NULL) {
        printf("Allocation of vectors failed!");
        return;
    }


    // Free matrices and vectors
    host_free_2d(Sigma_h);
    host_free_2d(L_h);
    //device_free_2d(Sigma_d, Sigma_log);
    //device_free_2d(L_d, L_log);
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(z_d);
    cudaFreeHost(z_h);




}