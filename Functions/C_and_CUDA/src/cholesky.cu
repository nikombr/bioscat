

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
extern "C" {
#include "../lib/GaussianProcess.h"

__device__ double reduction_step(double **L, int i, int k, int s, int idx) {
    
    double reg = idx < s ? L[i][idx]*L[k][idx] : 0; // Make sure all the threads have indeed a value - should result in less error
    for (int dist = 16; dist > 0; dist /= 2)
        reg += __shfl_down_sync(-1, reg, dist);

    return reg;

}

__device__ double reduction_step2(double *array, int n, int idx) {
    
    double reg = idx < n ? array[idx] : 0; // Make sure all the threads have indeed a value - should result in less error
    for (int dist = 16; dist > 0; dist /= 2)
        reg += __shfl_down_sync(-1, reg, dist);

    return reg;

}


__global__ void  sum(double **L, int i, int k, int s, double *res) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    double reg = reduction_step(L, i, k, s, idx);

    // Allocate shared memory statically
    __shared__ double smem[32];

    if (threadIdx.x % 32 == 0) smem[threadIdx.x/32] = reg;
    __syncthreads();

    idx = threadIdx.x;
    reg = reduction_step2(smem, 32, idx);

    if (threadIdx.x == 0) atomicAdd(res, reg);
    

}

__global__ void update_diagonal(double ** Sigma, double **L, int k, double ell2) {

    

}

__global__ void update(double ** Sigma, double **L, int k, int i, double ell2) {

    

}


void GaussianProcess::cholesky() {

    if (device) {

        for (int k = 0; k < n; k++) {
            for (int i = k + 1; i < n; i++) {

            }
        }

    }
    else {
        double tempsum;
        for (int k = 0; k < n; k++) {
            tempsum = 0;
            for (int s = 0; s < k-1; s++) {
                tempsum += L_h[k][s]*L_h[k][s];
            }
            L_h[k][k] = sqrt(Sigma_h[k][k] - tempsum);
            for (int i = k; i < n; i++) {
                tempsum = 0;
                for (int s = 0; s < k-1; s++) {
                    tempsum += L_h[i][s]*L_h[k][s];
                }
                L_h[k][k] = (Sigma_h[i][k] - tempsum)/L_h[k][k];

            }
        }

        
    }
    
}

}