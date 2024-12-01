
#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../../lib/BioScat.h"
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/combinePolarisation.h"

extern "C" {
using namespace std;

void BioScat::computeScatteredFields() {
    computeScatteredFields(beta);
}

void BioScat::computeScatteredFields(double beta) {
    combinePolarisation(E_scat_pol, E_scat, beta, deviceComputation);
    combinePolarisation(H_scat_pol, H_scat, beta, deviceComputation);
    
}

void BioScat::computeInteriorFields() {
    computeInteriorFields(beta);
}

void BioScat::computeInteriorFields(double beta) {
    combinePolarisation(E_int_pol, E_int, beta, deviceComputation);
    combinePolarisation(H_int_pol, H_int, beta, deviceComputation);
    
}

void BioScat::computeIncidentFields() {
    computeIncidentFields(beta);
}

void BioScat::computeIncidentFields(double beta) {

    combinePolarisation(E_inc_pol, E_inc, beta, deviceComputation);
    combinePolarisation(H_inc_pol, H_inc, beta, deviceComputation);
}

}