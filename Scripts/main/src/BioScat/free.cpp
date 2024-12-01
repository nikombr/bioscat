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


void BioScat::free() {

    for (int i = 0; i < num_segments; i++) {
        segments[i].free();
    }

    delete[] segments;
    reflectance.free();
    x_obs.free();
    y_obs.free();

    E_scat.free();
    H_scat.free();
    E_inc.free();
    H_inc.free();
    E_int.free();
    H_int.free();
    for (int i = 0; i < 2; i++) {
        E_scat_pol[i].free();
        H_scat_pol[i].free();
        E_inc_pol[i].free();
        H_inc_pol[i].free();
        E_int_pol[i].free();
        H_int_pol[i].free();
    }
    cudaDeviceReset();
    if (printOutput) printf("DESTRUCTED!\n");

}

}