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

void BioScat::reset() {
    if (deviceComputation) {
        // Initialize all arrays to zero on the host
        for (int i = 0; i < 2; i++) {
            E_scat_pol[i].setDeviceZero();
            H_scat_pol[i].setDeviceZero();
            E_inc_pol[i].setDeviceZero();
            H_inc_pol[i].setDeviceZero();
            E_int_pol[i].setDeviceZero();
            H_int_pol[i].setDeviceZero();
        }
    }
    else {
        // Initialize all arrays to zero on the host
        for (int i = 0; i < 2; i++) {
            E_scat_pol[i].setHostZero();
            H_scat_pol[i].setHostZero();
            E_inc_pol[i].setHostZero();
            H_inc_pol[i].setHostZero();
            E_int_pol[i].setHostZero();
            H_int_pol[i].setHostZero();
        }
    }
}

void BioScat::resetFarField() {
    if (deviceComputation) {
        // Initialize all arrays to zero on the host
        for (int i = 0; i < 2; i++) {
            F_pol[i].setDeviceZero();
        }
    }
    else {
        // Initialize all arrays to zero on the host
        for (int i = 0; i < 2; i++) {
            F_pol[i].setHostZero();
        }
    }
}
}