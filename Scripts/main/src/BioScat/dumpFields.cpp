#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../lib/BioScat.h"
#include "../lib/Segment.h"
#include "../lib/RealMatrix.h"
#include "../lib/combinePolarisation.h"
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
extern "C" {
using namespace std;

void BioScat::dumpFields() {

    if (deviceComputation) {
        E_scat.toHost();
        H_scat.toHost();
        E_inc.toHost();
        H_inc.toHost();
    }
    

    char filename[256];

    // Save scattered electric fields
    sprintf(filename, "../../../Results/forward/Ex_scat.txt");
    E_scat.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ey_scat.txt");
    E_scat.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ez_scat.txt");
    E_scat.z.dumpResult(filename);

    // Save scattered magnetic fields
    sprintf(filename, "../../../Results/forward/Hx_scat.txt");
    H_scat.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hy_scat.txt");
    H_scat.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hz_scat.txt");
    H_scat.z.dumpResult(filename);

    // Save incident electric fields
    sprintf(filename,"../../../Results/forward/Ex_inc.txt");
    E_inc.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ey_inc.txt");
    E_inc.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ez_inc.txt");
    E_inc.z.dumpResult(filename);

    // Save incident magnetic fields
    sprintf(filename,"../../../Results/forward/Hx_inc.txt");
    H_inc.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hy_inc.txt");
    H_inc.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hz_inc.txt");
    H_inc.z.dumpResult(filename);
}
}