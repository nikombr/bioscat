
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
extern "C" {
#include "../lib/GaussianProcess.h"

__device__ __host__ double squared_exponential(double* a, double* b, int dim,  double* hyper) {
    double tau = hyper[0];
    double ell = hyper[1];
    double temp = 0;
    for (int i = 0; i < dim; i++) {
        temp += (a[i]-b[i]) * (a[i]-b[i]);
    }

    return tau*tau*exp(-temp/(2*ell*ell));

}


#define kfunc squared_exponential

__global__ void covariance_matrix_device_1d(double ** Sigma, double* x, int dim, double* hyper, int n) {
    double a[1], b[1];
    double val;
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    if (k < n && i >= k && i < n) {
        a[0] = x[k];
        b[0] = x[i];
        val =  kfunc(a, b, dim, hyper);
        Sigma[i][k] = val;
        Sigma[k][i] = val;
        if (k == i) Sigma[i][k] += 1e-10;
    }
}

__global__ void covariance_matrix_device_2d(double ** Sigma, double* x, double* y, int dim, double* hyper, int n) {
    double val;
    double a[2], b[2];
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    int i = threadIdx.y + blockIdx.y * blockDim.y;
    if (k < n && i >= k && i < n) {
        a[0] = x[k];
        b[0] = x[i];
        a[1] = y[k];
        b[1] = y[i];
        val =  kfunc(a, b, dim, hyper);
        Sigma[i][k] = val;
        Sigma[k][i] = val;
        if (k == i) Sigma[i][k] += 1e-10;
    }
}



void GaussianProcess::covariance_matrix() {

    if (device) {
        if (dim == 1) {
            // Blocks and threads
            dim3 dimBlock(32,32);
            dim3 dimGrid((n+dimBlock.x-1)/dimBlock.x,(n+dimBlock.y-1)/dimBlock.y);
            // Call kernel
            covariance_matrix_device_1d<<<dimGrid, dimBlock>>>(M_d, x_d, dim, hyper_d, n);
            cudaDeviceSynchronize();
        }
        else if (dim == 2) {
            // Blocks and threads
            dim3 dimBlock(32,32);
            dim3 dimGrid((n+dimBlock.x-1)/dimBlock.x,(n+dimBlock.y-1)/dimBlock.y);
            // Call kernel
            covariance_matrix_device_2d<<<dimGrid, dimBlock>>>(M_d, x_d, y_d, dim, hyper_d, n);
            cudaDeviceSynchronize();

        }


    }
    else {
        double val;

        if (dim == 1) {
            double a[1], b[1];
            for (int k = 0; k < n; k++) {
                for (int i = k; i < n; i++) {
                    a[0] = x_h[k];
                    b[0] = x_h[i];
                    val = kfunc(a, b, dim, hyper_h);
                    M_h[i][k] = val;
                    M_h[k][i] = val;
                    if (k == i) M_h[i][k] += 1e-10;
                }
            }
        }
        else if (dim == 2) {
            double a[2], b[2];
            for (int k = 0; k < n; k++) {
                for (int i = k; i < n; i++) {
                    
                    a[0] = x_h[k];
                    b[0] = x_h[i];
                    a[1] = y_h[k];
                    b[1] = y_h[i];

                    val = kfunc(a, b, dim, hyper_h);
                    M_h[i][k] = val;
                    M_h[k][i] = val;
                    if (k == i) M_h[i][k] += 1e-10;
                }
            }
        }


    }

}

}