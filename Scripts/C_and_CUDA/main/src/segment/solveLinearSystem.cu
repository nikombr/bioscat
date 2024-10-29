#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;


// LAPACK routine for solving linear system
void dgels_(const char * trans, const int * m, const int * n, const int * nrhs, double * A, const int * lda, double * B,  const int * ldb, double * work, int * lwork,int * info);
//#include <lapacke.h>
void Segment::solveLinearSystem() {


    

    /*A = RealMatrix(2*(n_ext + n_int),4 * num_test_points);
    b = RealMatrix(4 * num_test_points);

    for (int r = 0; r < 2*num_test_points; r++) {
        for (int c = 0; c < n_ext + n_int; c++) {
            A.setHostValue(c,                   r,                     A_real.getHostValue(r,c));
            A.setHostValue(c + n_ext + n_int, r,                      - A_imag.getHostValue(r,c));
            A.setHostValue(c,         r + 2*num_test_points,           A_imag.getHostValue(r,c));
            A.setHostValue(c + n_ext + n_int,r + 2*num_test_points,    A_real.getHostValue(r,c));
        }
    }
    for (int r = 0; r < 2*num_test_points; r++) {
        b.setHostValue(r,                     b_real.getHostValue(r));
        b.setHostValue(r + 2*num_test_points, b_imag.getHostValue(r));
    }*/
    /*
    double Atest[12] = {1, 2,  3,  4,
                        5, 6,  7,  8,
                        9, 10, 11, 12};
    */
    /*
    double Atest[12] = {1,5,9,
                        2,6,10}; 
    */
    double Atest[12] = {1, 2, 
                        5, 6, 
                        9, 10}; 

    double btest[4];
    btest[0] = 5;
    btest[1] = 17;
    btest[2] = 29;
    int rows = 3;
    int cols = 2;

    char trans = 'T';
    int m = A.cols;//4 * num_test_points;
    int n = A.rows;//2*(n_ext + n_int);
    //printf("%d %d\n",m,n);
    int nrhs = 1; 
    int lda = m;
    int ldb = std::max(m, n);
    int info;
    double work_query;
    int lwork = -1;
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            //printf("%.6e\t",A.getHostValue(j,i));
        }
        //printf("\n");
    }


    dgels_(&trans, &m, &n, &nrhs, A.getHostPointer(), &lda, b.getHostPointer(), &ldb, &work_query, &lwork, &info);
    
    lwork = (int)work_query;
    double *work = (double*)malloc(lwork * sizeof(double));

    /*printf("\n\n");
    printf("A:\n");
    for (int j = 0; j < 2*num_test_points; j++) {
        printf("%f\t + i(%f)\n",A_real.getHostValue(j,0),A_imag.getHostValue(j,0));
    }*/

    printf("Solving linear system now.\n");
    dgels_(&trans, &m, &n, &nrhs, A.getHostPointer(), &lda, b.getHostPointer(), &ldb, work, &lwork, &info);
    //printf("HEJ\n");
    //info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, A, lda, B, ldb);
    if (info != 0) {
        printf("An error occurred in solving: %d\n", info);
    }

    C = ComplexMatrix(n_int);
    D = ComplexMatrix(n_ext);

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
    x_test_top.free();      y_test_top.free();
    x_test_right.free();    y_test_right.free();
    x_test_bottom.free();   y_test_bottom.free();
    x_test_left.free();     y_test_left.free();
    n_x.free();             n_y.free();

    /*printf("b:\n");

    for (int j = 0; j < 4; j++) {
        printf("%f\n",btest[j]);
    }

    printf("C:\n");

    for (int j = 0; j < n_int; j++) {
        printf("%e\t + i(%e)\n",b.getHostValue(j),b.getHostValue(j+n_ext+n_int));
    }

    printf("D:\n");

    for (int j = n_int; j < n_int + n_ext; j++) {
        printf("%e\t + i(%e)\n",b.getHostValue(j),b.getHostValue(j+n_ext+n_int));
    }*/
    
 
}


}