// Part of segment that applies both to 2D and 3D
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/2D/Segment.h"
#include "../lib/RealMatrix.h"
using namespace std;

void Segment::computeIncidentFieldVectors(RealMatrix y) {
    
    int rows = y.rows;

    bool Ex_bool = scenario == 1 ? false : true;
    bool Ez_bool = scenario == 1 ? true  : false;
    bool Hx_bool = scenario == 1 ? true  : false;
    bool Hz_bool = scenario == 1 ? false : true;
    bool Ey_bool = false;
    bool Hy_bool = false;

    E_inc_vector = Field(rows, Ex_bool, Ey_bool, Ez_bool);
    H_inc_vector = Field(rows, Hx_bool, Hy_bool, Hz_bool);

    if (scenario == 1) {
        for (int j = 0; j < rows; j++) E_inc_vector.z.setHostRealValue(j,                  cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_inc_vector.z.setHostImaginaryValue(j,               sin(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.x.setHostRealValue(j,    -1/constants.eta0*cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.x.setHostImaginaryValue(j, -1/constants.eta0*sin(constants.k0*y.getHostValue(j)));
    }
    else {
        for (int j = 0; j < rows; j++) E_inc_vector.x.setHostRealValue(j,                 cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_inc_vector.x.setHostImaginaryValue(j,              sin(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.z.setHostRealValue(j,    1/constants.eta0*cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.z.setHostImaginaryValue(j, 1/constants.eta0*sin(constants.k0*y.getHostValue(j)));
    }


}

void Segment::computeReflectedFieldVectors(RealMatrix y) {
    
    int rows = y.rows;

    bool Ex_bool = scenario == 1 ? false : true;
    bool Ez_bool = scenario == 1 ? true  : false;
    bool Hx_bool = scenario == 1 ? true  : false;
    bool Hz_bool = scenario == 1 ? false : true;
    bool Ey_bool = false;
    bool Hy_bool = false;

    E_ref_vector = Field(rows, Ex_bool, Ey_bool, Ez_bool);
    H_ref_vector = Field(rows, Hx_bool, Hy_bool, Hz_bool);

    if (scenario == 1) {
        for (int j = 0; j < rows; j++) E_ref_vector.z.setHostRealValue(j,                  constants.Gamma_ref*cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_ref_vector.z.setHostImaginaryValue(j,              -constants.Gamma_ref*sin(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_ref_vector.x.setHostRealValue(j,     1/constants.eta0*constants.Gamma_ref*cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_ref_vector.x.setHostImaginaryValue(j, -1/constants.eta0*constants.Gamma_ref*sin(constants.k0*y.getHostValue(j)));
    }
    else {
        for (int j = 0; j < rows; j++) E_ref_vector.x.setHostRealValue(j,                 constants.Gamma_ref*cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_ref_vector.x.setHostImaginaryValue(j,             -constants.Gamma_ref*sin(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_ref_vector.z.setHostRealValue(j,   -1/constants.eta0*constants.Gamma_ref*cos(constants.k0*y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_ref_vector.z.setHostImaginaryValue(j, 1/constants.eta0*constants.Gamma_ref*sin(constants.k0*y.getHostValue(j)));
    }


}
  

void Segment::setupRightHandSide() {
    // Setup right-hand side in linear system
    b = RealMatrix(4 * num_test_points);

    double val, shift;
    int imaginaryShift = 2*num_test_points;
    ComplexMatrix *firstField_inc, * firstField_ref, *secondField_inc,*secondField_ref;

    if (scenario == 1) {
        firstField_inc = &E_inc_vector.z;
        firstField_ref = &E_ref_vector.z;
        secondField_inc = &H_inc_vector.x;
        secondField_ref = &H_ref_vector.x;
    }
    else if (scenario == 2) {
        firstField_inc = &H_inc_vector.x;
        firstField_ref = &H_ref_vector.x;
        secondField_inc = &E_inc_vector.z;
        secondField_ref = &E_ref_vector.z;
    }
    else {
        printf("You have to choose either 1 or 2 as the scenario!\n");
        return;
    }

    for (int j = 0; j < num_test_points; j++) {
        val =  - firstField_inc->getHostRealValue(j) - firstField_ref->getHostRealValue(j);
        b.setHostValue(j, val);
        val =  - firstField_inc->getHostImaginaryValue(j) - firstField_ref->getHostImaginaryValue(j);
        b.setHostValue(j + imaginaryShift, val);
    }

    shift = num_test_points;
    for (int j = 0; j < n_top; j++) {
        val  = secondField_inc->getHostRealValue(j) + secondField_ref->getHostRealValue(j);
        val *= n_y.getHostValue(j);
        b.setHostValue(j + shift, val);
        val  = secondField_inc->getHostImaginaryValue(j) + secondField_ref->getHostImaginaryValue(j);
        val *= n_y.getHostValue(j);
        b.setHostValue(j + shift + imaginaryShift, val);
    }

    shift += n_top;
    for(int j = 0; j < n_right; j++) {
        val = 0.0;
        b.setHostValue(j + shift, val);
        b.setHostValue(j + shift + imaginaryShift, val);
    }

    shift += n_right;
    for(int j = 0; j < n_bottom; j++) {
        val  = secondField_inc->getHostRealValue(j) + secondField_ref->getHostRealValue(j);
        b.setHostValue(j + shift, val);
        val  = secondField_inc->getHostImaginaryValue(j) + secondField_ref->getHostImaginaryValue(j);
        b.setHostValue(j + shift + imaginaryShift, val);
    }

    shift += n_bottom;
    for(int j = 0; j < n_left; j++) {
        val = 0.0;
        b.setHostValue(j + shift, val);
        b.setHostValue(j + shift + imaginaryShift, val);
    }
    
    return;
 
}

// LAPACK routine for solving linear system
void dgels_(const char * trans, const int * m, const int * n, const int * nrhs, double * A, const int * lda, double * B,  const int * ldb, double * work, int * lwork,int * info);
//#include <lapacke.h>
void Segment::solveLinearSystem() {
    A = RealMatrix(4 * num_test_points, 2*(n_ext + n_int));

    for (int r = 0; r < 2*num_test_points; r++) {
        for (int c = 0; c < n_ext + n_int; c++) {
            A.setHostValue(r,                     c,                   A_real.getHostValue(r,c));
            A.setHostValue(r,                     c + n_ext + n_int, - A_imaginary.getHostValue(r,c));
            A.setHostValue(r + 2*num_test_points, c,                   A_imaginary.getHostValue(r,c));
            A.setHostValue(r + 2*num_test_points, c + n_ext + n_int,   A_real.getHostValue(r,c));
        }
    }

    char trans = 'T';
    int m = 4 * num_test_points;
    int n = 2*(n_ext + n_int);
    printf("%d %d\n",m,n);
    int nrhs = 1; 
    int lda = n;
    int ldb = m;
    int info;
    double work_query;
    int lwork = -1;
    dgels_(&trans, &n, &m, &nrhs, A.getHostPointer(), &lda, b.getHostPointer(), &ldb, &work_query, &lwork, &info);

    lwork = (int)work_query;
    double *work = (double*)malloc(lwork * sizeof(double));

    for (int j = 0; j < 4; j++) {
        printf("%e\n",b.getHostValue(j));
    }

    dgels_(&trans, &n, &m, &nrhs, A.getHostPointer(), &lda, b.getHostPointer(), &ldb, work, &lwork, &info);
    //info = LAPACKE_dgels(LAPACK_ROW_MAJOR, 'N', m, n, nrhs, A, lda, B, ldb);
    if (info != 0) {
        printf("An error occurred in solving: %d\n", info);
    }

    for (int j = 0; j < n; j++) {
        printf("%e\n",b.getHostValue(j));
    }
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            //printf("%e\n",A.getHostValue(j,i));
        }
    }
}


}