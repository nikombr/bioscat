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

void BioScat::prepareForward(double lambda) {
    for (int i = 0; i < num_segments; i++) {
        segments[i].newWavelength(lambda);
    }
}

void BioScat::prepareForward(double beta, double lambda) {
    this->beta = beta;
    for (int i = 0; i < num_segments; i++) {
        segments[i].newWavelength(lambda);
    }
}

}