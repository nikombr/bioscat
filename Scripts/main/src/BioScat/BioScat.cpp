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

BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points) {

    this->protein_structure = protein_structure;
    this->num_segments = num_segments;
    this->total_grid_points = total_grid_points;

}

BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points, bool deviceComputation) {

    this->protein_structure = protein_structure;
    this->num_segments = num_segments;
    this->total_grid_points = total_grid_points;
    this->deviceComputation = deviceComputation;

}

}