#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "../../lib/utils/ComplexMatrix.h"
#include "../../lib/BioScat.h"
extern "C" {
using namespace std;


void forward(double *x, double*y, int n ,char * protein_structure, int num_segments, int total_grid_points, double beta, double lambda, int deviceComputation_int, int printOutput_int) {

    bool deviceComputation = deviceComputation_int == 1 ? true : false;

    bool printOutput = printOutput_int == 1 ? true : false;

    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points, deviceComputation, x, y, n, printOutput);

    bioscat.getNanostructure();

    bioscat.setupSegments();
    bool compute_interior = false;
    if (bioscat.y_obs.findMin() < bioscat.nanostructure.f.findMax()) {
        compute_interior = true;
        bioscat.computeInterior = true;
    }
   

    bioscat.reset();

    bioscat.prepareForward(beta, lambda);

    for (int polarisation = 1; polarisation <= 2; polarisation++) {
       
        bioscat.forwardSolver(polarisation);
        
        bioscat.computeScatteredSubFields();
        if (compute_interior) bioscat.computeInteriorSubFields();
        //bioscat.computeReflectedSubFields();
        bioscat.computeIncidentSubFields();

    }

    bioscat.computeScatteredFields();
    if (compute_interior) bioscat.computeInteriorFields();

    //bioscat.computeReflectedFields();

    bioscat.computeIncidentFields();

    bioscat.dumpFields();

    bioscat.free();

}



}