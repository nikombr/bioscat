#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../../lib/BioScat.h"
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
extern "C" {
using namespace std;

void BioScat::setupObservationPoints(double *x, double*y, int n) {
    n_obs = n;
    x_obs = RealMatrix(n);
    y_obs = RealMatrix(n);
    for (int i = 0; i < n; i++) x_obs.setHostValue(i, x[i]);
    for (int i = 0; i < n; i++) y_obs.setHostValue(i, y[i]);
    x_obs.toDevice();
    y_obs.toDevice();

    // Allocate fields
    E_scat = Field(n);
    H_scat = Field(n);
    E_int = Field(n);
    H_int = Field(n);
    E_inc  = Field(n);
    H_inc  = Field(n);
    for (int i = 0; i < 2; i++) {
        E_scat_pol[i] = Field(n);
        H_scat_pol[i] = Field(n);
        E_int_pol[i]  = Field(n);
        H_int_pol[i]  = Field(n);
        E_inc_pol[i]  = Field(n);
        H_inc_pol[i]  = Field(n);
    }

    // Allocate reflectance array
    //reflectance = RealMatrix(n);
}

void BioScat::setupObservationPoints(double *phi, int n) {
    n_obs = n;
    phi_obs = RealMatrix(n);
    for (int i = 0; i < n; i++) phi_obs.setHostValue(i, phi[i]);
    phi_obs.toDevice();

    // Allocate fields
    F = RealMatrix(n);
    for (int i = 0; i < 2; i++) {
        F_pol[i] = ComplexMatrix(n);
    }

    // Allocate reflectance array
    //reflectance = RealMatrix(n);
}

}