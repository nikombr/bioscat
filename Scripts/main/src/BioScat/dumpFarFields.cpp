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

void BioScat::dumpFarFields() {

    if (deviceComputation) {
        F.toHost();
    }

    char filename[256];
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");

    // Save scattered electric fields
    sprintf(filename, "%s/Results/forward/farFieldPattern.txt", dir);
    F.dumpResult(filename);
}
}