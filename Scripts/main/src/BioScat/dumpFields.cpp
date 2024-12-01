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

void BioScat::dumpFields() {

    if (deviceComputation) {
        E_scat.toHost();
        H_scat.toHost();
        E_inc.toHost();
        H_inc.toHost();
        E_int.toHost();
        H_int.toHost();
    }
    

    char filename[256];

    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");

    // Save scattered electric fields
    sprintf(filename, "%s/Results/forward/Ex_scat.txt", dir);
    E_scat.x.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Ey_scat.txt", dir);
    E_scat.y.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Ez_scat.txt", dir);
    E_scat.z.dumpResult(filename);

    // Save scattered magnetic fields
    sprintf(filename, "%s/Results/forward/Hx_scat.txt", dir);
    H_scat.x.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Hy_scat.txt", dir);
    H_scat.y.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Hz_scat.txt", dir);
    H_scat.z.dumpResult(filename);

    if (computeInterior) {
        // Save interior electric fields
        sprintf(filename, "%s/Results/forward/Ex_int.txt", dir);
        E_int.x.dumpResult(filename);
        sprintf(filename,"%s/Results/forward/Ey_int.txt", dir);
        E_int.y.dumpResult(filename);
        sprintf(filename,"%s/Results/forward/Ez_int.txt", dir);
        E_int.z.dumpResult(filename);

        // Save interior magnetic fields
        sprintf(filename, "%s/Results/forward/Hx_int.txt", dir);
        H_int.x.dumpResult(filename);
        sprintf(filename,"%s/Results/forward/Hy_int.txt", dir);
        H_int.y.dumpResult(filename);
        sprintf(filename,"%s/Results/forward/Hz_int.txt", dir);
        H_int.z.dumpResult(filename);
    }

    // Save incident electric fields
    sprintf(filename,"%s/Results/forward/Ex_inc.txt", dir);
    E_inc.x.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Ey_inc.txt", dir);
    E_inc.y.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Ez_inc.txt", dir);
    E_inc.z.dumpResult(filename);

    // Save incident magnetic fields
    sprintf(filename,"%s/Results/forward/Hx_inc.txt", dir);
    H_inc.x.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Hy_inc.txt", dir);
    H_inc.y.dumpResult(filename);
    sprintf(filename,"%s/Results/forward/Hz_inc.txt", dir);
    H_inc.z.dumpResult(filename);
}
}