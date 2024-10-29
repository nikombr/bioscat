#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
#include "../lib/BioScat.h"
using namespace std;


void forward(double *x, double*y, int n ,char * protein_structure, int num_segments) {

    BioScat bioscat = BioScat(protein_structure, num_segments);

    bioscat.getNanostructure();

    bioscat.getSegments();

    bioscat.setupObservationPoints(x, y, n);

    bioscat.forwardSolver(1);

    bioscat.computeScatteredFields();

    bioscat.computeReflectedFields();

    bioscat.computeIncidentFields();

    printf("Hej fra C\n");

}



}