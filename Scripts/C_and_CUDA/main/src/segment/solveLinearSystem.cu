#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <cublas_v2.h>
#include <omp.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

// LAPACK routine for solving linear system
void dgels_(const char * trans, const int * m, const int * n, const int * nrhs, double * A, const int * lda, double * B,  const int * ldb, double * work, int * lwork,int * info);

// Need to transpose manually as cublas dgels does not seem to have implemented this yet
__global__ void transpose(double * input, double * output, int rows, int cols) {
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    if (r < rows && c < cols) {
        output[c*rows + r] = input[r*cols + c];
    }
}

void Segment::solveLinearSystem() {

    //C = ComplexMatrix(n_int);
    //D = ComplexMatrix(n_ext);

    char trans;
    int m, n, nrhs, lda, ldb, info, lwork;
    double work_query;
    double *work;

    /*// Test if it works
    double Atest_h[12] = {1, 2, 3,
                        5, 7, 7,
                        9, 10, 11,
                        12, 13, 14}; 
    //double Atest_h[12] = {1, 5, 9, 12, 2, 7, 10, 13, 3, 7, 11, 14}; 

    double btest_h[4];
    btest_h[0] = 14;
    btest_h[1] = 40;
    btest_h[2] = 62;
    btest_h[3] = 80;
    int rows = 4;
    int cols = 3;
    
    m = rows;
    n = cols;
    nrhs = 1; 
    lda = m;
    ldb = std::max(m, n);
    printf("(cols, rows, ldb) = (%d, %d, %d)\n", cols, rows, ldb);
    lwork = -1;

    // cuBLAS handle creation
    cublasHandle_t handle;
    cublasStatus_t status;
    status = cublasCreate(&handle);

    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("cuBLAS initialization failed %d\n",status);
        return;
    }
    double * Atest_d, *Atest_T_d, *btest_d;

    cudaMalloc((void **) &Atest_d,     rows * cols * sizeof(double));
    cudaMalloc((void **) &Atest_T_d,     rows * cols * sizeof(double));
    cudaMalloc((void **) &btest_d,     rows * sizeof(double));
    cudaMemcpy(Atest_d,    Atest_h,    rows * cols * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(btest_d,    btest_h,    rows * sizeof(double), cudaMemcpyHostToDevice);

    // Blocks and threads
    dim3 dimBlock(32,32);
    dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x, (cols + dimBlock.y - 1)/dimBlock.y);
    transpose<<< dimGrid, dimBlock>>>(Atest_d, Atest_T_d, rows, cols);

    double * A_ptr_h = Atest_T_d; // have tried Atest_h
    double * b_ptr_h = btest_d;
    double ** A_ptr_d;
    double ** b_ptr_d;
    cudaMalloc(&A_ptr_d, sizeof(double*));
    cudaMalloc(&b_ptr_d, sizeof(double*));
    cudaMemcpy(A_ptr_d, &A_ptr_h, sizeof(double*), cudaMemcpyHostToDevice);
    cudaMemcpy(b_ptr_d, &b_ptr_h, sizeof(double*), cudaMemcpyHostToDevice);

    int * info_h;
    int* info_d;
    //cudaMalloc((void**)&info_d, sizeof(int));
    printf("HEJ FRA HER!\n");
    //int *devInfoArray; 
    //cudaMalloc((void**)&devInfoArray,  sizeof(int));
    cudaMallocHost((void **) &info_h,  sizeof(int));
    cudaMalloc((void **) &info_d,     sizeof(int));

    int devInfoArray[1] = { 0 };
        cudaDeviceSynchronize();
    status = cublasDgelsBatched(handle, CUBLAS_OP_N, m, n, nrhs, A_ptr_d, lda, b_ptr_d, ldb, info_h, NULL, 1);
    cudaMemcpy(btest_h,    btest_d,    rows * sizeof(double), cudaMemcpyDeviceToHost);

    //cudaMemcpy(info_h, info_d, sizeof(int), cudaMemcpyDeviceToHost);
    printf("HEJ FRA HER 2!\n");
    cudaDeviceSynchronize();

    printf("b:\n");

    for (int j = 0; j < 4; j++) {
        printf("%f\n",btest_h[j]);
    }

    // Check if cublasDtrmv was successful
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("cublasDgelsBatched failed with error code: %d\n", status);
        printf("info = %d\n",*info_h);
    }
    printf("info = %d\n",*info_h);

    // Destroy cuBLAS handle
    status = cublasDestroy(handle);
    if (status != CUBLAS_STATUS_SUCCESS) {
        printf("cuBLAS destruction failed\n");
        return;
    }
    

    /*dgels_(&trans, &m, &n, &nrhs, Atest, &lda, btest, &ldb, &work_query, &lwork, &info);
    
    lwork = (int)work_query;
    work = (double*)malloc(lwork * sizeof(double));

    printf("Solving linear system now.\n");
    dgels_(&trans, &m, &n, &nrhs, Atest, &lda, btest, &ldb, work, &lwork, &info);

    printf("b:\n");

    for (int j = 0; j < 4; j++) {
        printf("%f\n",btest[j]);
    }
    double Atest2[12] = {1, 2, 3, 5, 6, 7, 9, 10, 11,12, 13, 14}; 
    printf("Ab:\n");
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            printf("%f ",Atest2[r*cols+c]);
        }
        printf("\n");
    }
    printf("Ab:\n");
    for (int r = 0; r < rows; r++) {
        double val = 0.0;
        for (int c = 0; c < cols; c++) {
            printf("%d ",r*cols+c);
            val += Atest2[r*cols+c]*btest[c];
            printf("%f ",val);
        }
        printf("%f\n",val);
    }
/*printf("b:\n");
    for (int j = 0; j < n_int; j++) {
        printf("%e\n",b.getHostValue(j));
    }*/

    trans = 'T';
    m = A.cols;
    n = A.rows;
    nrhs = 1; 
    lda = m;
    ldb = std::max(m, n);
    lwork = -1;

    if (false) {
        // cuBLAS handle creation
        cublasHandle_t handle;
        cublasStatus_t status;
        status = cublasCreate(&handle);

        if (status != CUBLAS_STATUS_SUCCESS) {
            printf("cuBLAS initialization failed %d\n",status);
            return;
        }

        A.toDevice();
        b.toDevice();

        m = A.rows;
        n = A.cols;
         nrhs = 1; 
    lda = m;
    ldb = std::max(m, n);
    lwork = -1;
        double * A_T_d;
        cudaMalloc((void **) &A_T_d, A.rows * A.cols * sizeof(double));

        // Blocks and threads
        dim3 dimBlock(32,32);
        dim3 dimGrid((A.rows + dimBlock.x - 1)/dimBlock.x, (A.cols + dimBlock.y - 1)/dimBlock.y);
        transpose<<< dimGrid, dimBlock>>>(A.getDevicePointer(), A_T_d, A.rows, A.cols);

    
        double * A_ptr_h = A_T_d;
        double * b_ptr_h = b.getDevicePointer();
        double ** A_ptr_d;
        double ** b_ptr_d;
        cudaMalloc(&A_ptr_d, sizeof(double*));
        cudaMalloc(&b_ptr_d, sizeof(double*));
        cudaMemcpy(A_ptr_d, &A_ptr_h, sizeof(double*), cudaMemcpyHostToDevice);
        cudaMemcpy(b_ptr_d, &b_ptr_h, sizeof(double*), cudaMemcpyHostToDevice);


        int * info_h;
    int* info_d;
    //cudaMalloc((void**)&info_d, sizeof(int));
    printf("HEJ FRA HER!\n");
    //int *devInfoArray; 
    //cudaMalloc((void**)&devInfoArray,  sizeof(int));
    cudaMallocHost((void **) &info_h,  sizeof(int));
    cudaMalloc((void **) &info_d,     sizeof(int));
    double start = omp_get_wtime();
        cudaDeviceSynchronize();
        
        status = cublasDgelsBatched(handle, CUBLAS_OP_N, m, n, nrhs, A_ptr_d, lda, b_ptr_d, ldb, info_h, NULL, 1);
        
        cudaDeviceSynchronize();
        double end = omp_get_wtime();
        printf("time = %f\n",end-start);
        printf("HEJ FRA HER 2!\n");
        // Check if cublasDtrmv was successful
        if (status != CUBLAS_STATUS_SUCCESS) {
            printf("cublasDgelsBatched failed with error code: %d\n", status);
        }

        // Destroy cuBLAS handle
        status = cublasDestroy(handle);
        if (status != CUBLAS_STATUS_SUCCESS) {
            printf("cuBLAS destruction failed\n");
            return;
        }
        b.toHost();
    }
    else {
        dgels_(&trans, &m, &n, &nrhs, A.getHostPointer(), &lda, b.getHostPointer(), &ldb, &work_query, &lwork, &info);
    
        lwork = (int)work_query;
        work = (double*)malloc(lwork * sizeof(double));

        printf("Solving linear system now.\n");
        dgels_(&trans, &m, &n, &nrhs, A.getHostPointer(), &lda, b.getHostPointer(), &ldb, work, &lwork, &info);
        if (info != 0) {
            printf("An error occurred in solving: %d\n", info);
        }


    }
    for (int i = 0; i < n_int; i++) {
        C.setHostRealValue(i,b.getHostValue(i));
        C.setHostImagValue(i,b.getHostValue(i + n_ext + n_int));
    }
    
    for (int i = 0; i < n_ext; i++) {
        D.setHostRealValue(i,b.getHostValue(i + n_int));
        D.setHostImagValue(i,b.getHostValue(i + n_ext + 2*n_int));
    }


    // Free arrays that we no longer need
    //A.free();
    //b.free();
    //n_x.free();
    //n_y.free();
    //x_test.free();
    //y_test.free();

    /*printf("b:\n");

    for (int j = 0; j < 4; j++) {
        printf("%f\n",btest[j]);
    }*/

    /*printf("C:\n");

    for (int j = 0; j < n_int; j++) {
        printf("%e\t + i(%e)\n",C.getHostRealValue(j),C.getHostImagValue(j));
    }

    printf("D:\n");

    for (int j = 0; j < n_ext; j++) {
        printf("%e\t + i(%e)\n",D.getHostRealValue(j),D.getHostImagValue(j));
    }*/
    
 
}


}