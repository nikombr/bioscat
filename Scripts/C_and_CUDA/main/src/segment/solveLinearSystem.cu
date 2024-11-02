#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <cublas_v2.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

// LAPACK routine for solving linear system
void dgels_(const char * trans, const int * m, const int * n, const int * nrhs, double * A, const int * lda, double * B,  const int * ldb, double * work, int * lwork,int * info);

void Segment::solveLinearSystem() {

    //C = ComplexMatrix(n_int);
    //D = ComplexMatrix(n_ext);

    char trans;
    int m, n, nrhs, lda, ldb, info, lwork;
    double work_query;
    double *work;

    /*// Test if it works
    double Atest[12] = {1, 2, 3,
                        5, 7, 7,
                        9, 10, 11,
                        12, 13, 14}; 

    double btest[4];
    btest[0] = 14;
    btest[1] = 40;
    btest[2] = 62;
    btest[3] = 80;
    int rows = 4;
    int cols = 3;
    trans = 'T';
    m = cols;
    n = rows;
    nrhs = 1; 
    lda = m;
    ldb = std::max(m, n);
    lwork = -1;

    dgels_(&trans, &m, &n, &nrhs, Atest, &lda, btest, &ldb, &work_query, &lwork, &info);
    
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
    }*/
printf("b:\n");
    for (int j = 0; j < n_int; j++) {
        printf("%e\n",b.getHostValue(j));
    }

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

        double * A_ptr_h = A.getDevicePointer();
        double * b_ptr_h = b.getDevicePointer();
        double ** A_ptr_d;
        double ** b_ptr_d;
        cudaMalloc(&A_ptr_d, sizeof(double*));
        cudaMalloc(&b_ptr_d, sizeof(double*));
        cudaMemcpy(A_ptr_d, &A_ptr_h, sizeof(double*), cudaMemcpyHostToDevice);
        cudaMemcpy(b_ptr_d, &b_ptr_h, sizeof(double*), cudaMemcpyHostToDevice);

        //int info[1];
        int* info_d;
        cudaMalloc((void**)&info_d, sizeof(int));
        printf("HEJ FRA HER!\n");
        int *devInfoArray; 
        cudaMalloc((void**)&devInfoArray,  sizeof(int));
        
        status = cublasDgelsBatched(handle, CUBLAS_OP_T, m, n, nrhs, A_ptr_d, lda, b_ptr_d, ldb, &info, NULL, 1);
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
    A.free();
    b.free();
    n_x.free();
    n_y.free();
    x_test.free();
    y_test.free();

    /*printf("b:\n");

    for (int j = 0; j < 4; j++) {
        printf("%f\n",btest[j]);
    }*/

    printf("C:\n");

    for (int j = 0; j < n_int; j++) {
        printf("%e\t + i(%e)\n",C.getHostRealValue(j),C.getHostImagValue(j));
    }

    printf("D:\n");

    for (int j = 0; j < n_ext; j++) {
        printf("%e\t + i(%e)\n",D.getHostRealValue(j),D.getHostImagValue(j));
    }
    
 
}


}