#include <stdlib.h>
#include <cublas_v2.h>
#include <stdio.h>

extern "C" {
#include "../lib/GaussianProcess.h"
#include <cblas.h>
#include <ctime>

void GaussianProcess::generate_random_vector() {

    for (int i = 0; i < n; i++) {
        p_h[i] = ((double) rand())/((double) RAND_MAX);
    }

    if (n < 21) {
        printf("\n");
        for (int k = 0; k < n; k++) {
            if (k != n-1) printf("%.8f, ",p_h[k]);
            else printf("%.8f",p_h[k]);
            
        }
        printf("\n");
        printf("\n");
    }

    if (device) {
        // Send to device
        cudaMemcpy(p_d, p_h, n * sizeof(double), cudaMemcpyHostToDevice);

    }


    /*if (device) {

        // Blocks and threads
        dim3 dimBlock(32);
        dim3 dimGrid((n+dimBlock.x-1)/dimBlock.x);

        // Allocate memory for random number states
        curandState *states;
        cudaMalloc((void **)&states, n * sizeof(curandState));

        // Generate random states
        generate_random_states<<<dimGrid, dimBlock>>>(states, time(0), n);
        cudaDeviceSynchronize();

        // Generate random numbers
        generate_random_vector_device<<<dimGrid, dimBlock>>>(states, p_d, n);
        cudaDeviceSynchronize();


    }

    else {

        for (int i = 0; i < n; i++) {
            p_h[i] = ((double) rand())/((double) RAND_MAX);
        }

    }*/

}

void GaussianProcess::realisation() {

    // Seed the random number generator with the current time
    //srand(time(NULL));
    //srand(0);

    generate_random_vector();

    if (device) {

        // cuBLAS handle creation
        cublasHandle_t handle;
        cublasCreate(&handle);

        cublasDtrmv(handle, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_T, CUBLAS_DIAG_NON_UNIT, n, M_log, n, p_d, 1);

    }
    else {
        cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, *M_h, n, p_h, 1);
    }

}

// 0.84018772, 0.39438293, 0.78309922, 0.79844003, 0.91164736, 0.19755137

// 0.84018772, 0.82315774, 0.95280335, 0.90480132, 1.30937126, 0.97653560



}