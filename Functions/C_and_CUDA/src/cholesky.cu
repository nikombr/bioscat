

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cusolverDn.h>
extern "C" {
#include "../lib/GaussianProcess.h"
//#include <cblas.h>
//#include <cuda_runtime.h>
// LAPACK routine for Cholesky factorization
void dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);

/*__device__ double reduction_step(double **L, int i, int k, int idx) {
    
    double reg = idx < k ? L[i][idx]*L[k][idx] : 0; // Make sure all the threads have indeed a value - should result in less error
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


__global__ void  sum(double **L, int i, int k, double *res) {
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    double reg = reduction_step(L, i, k, idx);

    // Allocate shared memory statically
    __shared__ double smem[32];

    if (threadIdx.x % 32 == 0) smem[threadIdx.x/32] = reg;
    __syncthreads();

    idx = threadIdx.x;
    reg = reduction_step2(smem, 32, idx);

    if (threadIdx.x == 0) atomicAdd(res, reg);
    

}

__global__ void update_diagonal(double ** Sigma, double **L, int k, double *tempsum) {

    L[k][k] = sqrt(Sigma[k][k] - *tempsum);

}

__global__ void update(double ** Sigma, double **L, int k, int i, double *tempsum) {

    L[i][k] = (Sigma[k][i] - *tempsum)/L[k][k];

}

__global__ void init_zero(double *res) {
    *res = 0.0;
}*/

void GaussianProcess::cholesky() {

    if (device) {

        // Create cuSOLVER handle
        cusolverDnHandle_t cusolverH;
        cusolverDnCreate(&cusolverH);

        // Workspace size query
        int workspace_size = 0;
        cusolverDnDpotrf_bufferSize(cusolverH, CUBLAS_FILL_MODE_UPPER, n, M_log, n, &workspace_size);

        // Allocate workspace
        double *workspace;
        cudaMalloc((void**)&workspace, workspace_size * sizeof(double));

        // Allocate info variable
        int *info_d;
        cudaMalloc((void**)&info_d, sizeof(int));

        // Perform Cholesky factorization
        cusolverDnDpotrf(cusolverH, CUBLAS_FILL_MODE_UPPER, n, M_log, n, workspace, workspace_size, info_d);


        /*double *tempsum;
        cudaMalloc((void**)&tempsum, sizeof(double));
        init_zero<<<1, 1>>>(tempsum);
        cudaDeviceSynchronize();

        // Blocks and threads
        dim3 dimBlock(32);
        dim3 dimGrid((n+dimBlock.x-1)/dimBlock.x);

        for (int k = 0; k < n; k++) {

            init_zero<<<1, 1>>>(tempsum);
            cudaDeviceSynchronize();

            dimGrid.x = (k+dimBlock.x-1)/dimBlock.x;
            sum<<<dimGrid, dimBlock>>>(L_d, k, k, tempsum);
            cudaDeviceSynchronize();

            update_diagonal<<<1, 1>>>(Sigma_d, L_d, k, tempsum);
            cudaDeviceSynchronize();



            for (int i = k + 1; i < n; i++) {

                init_zero<<<1, 1>>>(tempsum);
                cudaDeviceSynchronize();

                sum<<<dimGrid, dimBlock>>>(L_d, i, k, tempsum);
                cudaDeviceSynchronize();

                update<<<1, 1>>>(Sigma_d, L_d, k, i, tempsum);
                cudaDeviceSynchronize();

            }
        }*/

    }
    else {
        char uplo = 'U';
        int N = n*n;
        int info;
        dpotrf_(&uplo, &n, *M_h, &n,&info);
        /*double tempsum;
        for (int k = 0; k < n; k++) {


            tempsum = 0;
            for (int s = 0; s < k; s++) {
                
                tempsum += L_h[k][s]*L_h[k][s];
            }

            L_h[k][k] = sqrt(Sigma_h[k][k] - tempsum);


            for (int i = k + 1; i < n; i++) {


                tempsum = 0;
                for (int s = 0; s < k; s++) {
                    tempsum += L_h[i][s]*L_h[k][s];
                }
                L_h[i][k] = (Sigma_h[k][i] - tempsum)/L_h[k][k];

            }
        }*/

        
    }
    
}

}