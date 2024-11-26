
#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../lib/Segment.h"
#include "../lib/BioScat.h"
#include "../lib/RealMatrix.h"
extern "C" {
using namespace std;


void BioScat::computeScatteredFields() {
    combinePolarisation(E_scat_pol, E_scat, beta, deviceComputation);
    combinePolarisation(H_scat_pol, H_scat, beta, deviceComputation);
    
}

void BioScat::computeIncidentFields() {

    combinePolarisation(E_inc_pol, E_inc, beta, deviceComputation);
    combinePolarisation(H_inc_pol, H_inc, beta, deviceComputation);
}

}