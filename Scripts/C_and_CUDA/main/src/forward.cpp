#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
#include "../lib/BioScat.h"
using namespace std;


void forward(char * protein_structure, int num_segments) {

    BioScat bioscat = BioScat(protein_structure, num_segments);

    bioscat.getNanostructure();

    bioscat.getSegments();

    bioscat.forwardSolver(1);

    printf("Hej fra C\n");

}



}