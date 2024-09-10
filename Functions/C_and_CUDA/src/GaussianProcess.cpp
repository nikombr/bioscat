
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
extern "C" {
#include "../lib/GaussianProcess.h"
#include "../lib/cudaMalloc2d.h"
using namespace std;

GaussianProcess::GaussianProcess(double* x, double* y, int n, double* hyper, int num, int dim) {

    if (dim == 1) {
        y = NULL;
    }

    // Save input parameters
    this->x_h       = x;
    this->y_h       = y;
    this->n         = n;
    this->hyper_h   = hyper;
    this->num       = num;
    this->dim       = dim;

    // Check if device is available
    int temp;
    cudaError_t cudaSuccess =  cudaGetDeviceCount(&temp);
    device = (!cudaSuccess && temp > 0) ? true : false;
    printf("devices = %d\n",cudaSuccess);

    string location = device ? "device" : "host";
    //device = false;

    if (y_h == NULL) printf("We are computing curves on %s!\n",location.c_str());
    else printf("We are computing planes on %s!\n",location.c_str());

    host_malloc_2d(&Sigma_h, n);
    printf("hej check\n");
    if (!device) {
        // Allocate matrices on host
        //host_malloc_2d(&Sigma_h, n);
        host_malloc_2d(&L_h, n);
        // Check allocation
        if (Sigma_h == NULL || L_h == NULL) {
            printf("Allocation of matrices failed on host!\n");
            return;
        }
        // Allocate vectors on host
        cudaMallocHost((void **) &p_h, n*sizeof(double));
        // Check allocation
        if (p_h == NULL) {
            printf("Allocation of vectors failed on host!\n");
            return;
        }
    }
    

    // Allocate vectors on host
    cudaMallocHost((void **) &z_h, n*sizeof(double));

    // Check allocation
    if (z_h == NULL) {
        printf("Allocation of vectors failed on host!\n");
        return;
    }
    
    if (device) {

        // Allocate matrices on device
        device_malloc_2d(&Sigma_d, &Sigma_log, n);
        device_malloc_2d(&L_d, &L_log, n);

        // Check allocation
        if (Sigma_d == NULL || L_d == NULL || Sigma_log == NULL || L_log == NULL) {
            printf("Allocation of matrices failed on device! %d\n",n);
            return;
        }

        // Allocate vectors on device
        cudaMalloc((void **) &x_d,     n*sizeof(double));
        if (y_h != NULL) cudaMalloc((void **) &y_d,     n*sizeof(double));
        cudaMalloc((void **) &z_d,     n*sizeof(double));
        cudaMalloc((void **) &p_d,     n*sizeof(double));
        cudaMalloc((void **) &hyper_d,     dim*sizeof(double));

        // Check allocation
        if (x_d == NULL || !(y_h == NULL || (y_h != NULL && y_d != NULL)) || z_d == NULL || p_d == NULL) {
            printf("Allocation of vectors failed on device!\n");
            return;
        }

        // Send to device
        cudaMemcpy(x_d, x_h, n * sizeof(double), cudaMemcpyHostToDevice);
        if (y_h != NULL) cudaMemcpy(y_d, y_h, n * sizeof(double), cudaMemcpyHostToDevice);

        // Send data
        cudaMemcpy(x_d, x_h, n * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(y_d, y_h, n * sizeof(double), cudaMemcpyHostToDevice);
        cudaMemcpy(hyper_d, hyper_h, dim * sizeof(double), cudaMemcpyHostToDevice);
        
    }
}

GaussianProcess::~GaussianProcess() {

    printf("start destruction!\n");

    if (device) {
        cudaMemcpy(*Sigma_h, Sigma_log, n * n * sizeof(double), cudaMemcpyDeviceToHost);
    }
    if (n < 21) {
        printf("print! %f\n", Sigma_h[0][0]);
        for (int k = 0; k < n; k++) {
            for (int i = 0; i < n; i++) {
                printf("%.4f ",Sigma_h[k][i]);
            }
            printf("\n");
        }
    }
    

    if (!device) {
        host_free_2d(Sigma_h);
        host_free_2d(L_h);
    }
    printf("Hej!\n");
    cudaFreeHost(x_h);
    if (dim == 2) cudaFreeHost(y_h);
    cudaFreeHost(z_h);
    cudaFreeHost(p_h);

    if (device) {
        device_free_2d(Sigma_d,Sigma_log);
        device_free_2d(L_d,L_log);
        cudaFree(x_d);
        if (dim == 2) cudaFree(y_d);
        cudaFree(z_d);
        cudaFree(p_d);
    }

    printf("Destructed!\n");

}

}