#include <stdlib.h>
#include <stdio.h>
#include "../lib/gaussian_process_inner.h"
#include <cuda_runtime_api.h>

void gaussian_process(double * x, double * y, int n, double * hyper, int num, int dim) {

    gaussian_process_inner(x, y, n, hyper, num, dim);

    /*// Initialize
    double ** Sigma_h, **L_h, **Sigma_d, **L_d;
    double *Sigma_log, *L_log, *z_h, *z_d, *x_d, *y_d, *p_d;

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
    cudaMalloc((void **) &p_d,     n*sizeof(double));
    cudaMallocHost((void **) &z_h, n*sizeof(double));

    // Check allocation
    if (x_d == NULL || y_d == NULL || x_d == NULL || z_h == NULL) {
        printf("Allocation of vectors failed!");
        return;
    }

    // Send to device
    cudaMemcpy(x_d, x_h, n * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(y_d, y_h, n * sizeof(double), cudaMemcpyHostToDevice);

    // Form covariance matrix

    // Form random vector

    // Compute matrix-vector-product

    // Send result to host - måske ikke nødvendigt, hvis vi bare vil skrive det til en fil og det kun afhænger af lokationen  
    
    cudaMemcpy(z_h, z_d, n * sizeof(double), cudaMemcpyDeviceToHost);

    // Free matrices and vectors
    host_free_2d(Sigma_h);
    host_free_2d(L_h);
    device_free_2d(Sigma_d, Sigma_log);
    device_free_2d(L_d, L_log);
    cudaFree(x_d);
    cudaFree(y_d);
    cudaFree(z_d);
    cudaFreeHost(z_h);*/

    printf("Hej fra C\n");



}