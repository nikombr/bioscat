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

void BioScat::setupSegments() {

    if (deviceComputation) {
        if (printOutput) printf("We are computing on the device!\n");
    }
    else {
        if (printOutput) printf("We are computing on the host!\n");
    }

    double start = omp_get_wtime();

    for (int i = 0; i < num_segments; i++) {
        segments[i].current_segment = i;
        segments[i].n_obs = n_obs;
        segments[i].setup(this->nanostructure, total_grid_points, num_segments);
    }


    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %e seconds to setup the segments!\n\n",end-start);

}

void BioScat::setupSegments(Nanostructure nanostructure) {

    double start = omp_get_wtime();

    for (int i = 0; i < num_segments; i++) {
        //segments[i] = Segment();
        segments[i].current_segment = i;
        segments[i].setup(nanostructure, total_grid_points, num_segments);
    }

    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %e seconds to setup the segments!\n\n",end-start);
}

}