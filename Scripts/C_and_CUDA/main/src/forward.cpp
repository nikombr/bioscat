#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
#include "../lib/BioScat.h"
using namespace std;


void forward(double *x, double*y, int n ,char * protein_structure, int num_segments, int total_grid_points, double beta, double lambda) {

    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points);
    bioscat.printOutput = true;

    bioscat.setupObservationPoints(x, y, n);

    bioscat.getNanostructure();

    bioscat.getSegments();

    bioscat.reset();

    bioscat.prepareForward(beta, lambda);

    for (int polarisation = 1; polarisation <= 2; polarisation++) {
       
        bioscat.forwardSolver(polarisation);
        
        bioscat.computeScatteredSubFields();
        bioscat.computeReflectedSubFields();
        bioscat.computeIncidentSubFields();

    }

    bioscat.computeScatteredFields();

    bioscat.computeReflectedFields();

    bioscat.computeIncidentFields();

    bioscat.dumpFields();

    bioscat.free();

}



}