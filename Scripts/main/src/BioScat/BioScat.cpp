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
#include "../../lib/combinePolarisation.h"
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
extern "C" {
using namespace std;

BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points) {

    this->protein_structure = protein_structure;
    this->num_segments = num_segments;
    this->total_grid_points = total_grid_points;

}

BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points, bool deviceComputation,double *phi_obs, int n_obs, bool printOutput) {

    this->protein_structure = protein_structure;
    this->num_segments      = num_segments;
    this->total_grid_points = total_grid_points;
    this->deviceComputation = deviceComputation;
    this->printOutput       = printOutput;

    this->n_obs = n_obs;
    this->phi_obs = RealMatrix(n_obs);
    for (int i = 0; i < n_obs; i++) this->phi_obs.setHostValue(i, phi_obs[i]);
    this->phi_obs.toDevice();

    // Allocate fields
    F = RealMatrix(n_obs);
    for (int i = 0; i < 2; i++) {
        F_pol[i] = ComplexMatrix(n_obs);
    }

    int sideSteps      = 0.75*total_grid_points;
    int segment_length = total_grid_points / num_segments;
    int n_top          = segment_length - 2;
    int n_bottom       = segment_length - 2;
    int n_right        = sideSteps - 1;
    int n_left         = sideSteps - 1;
    int n_int          = n_top + n_right + n_bottom + n_left - 16;
    int n_ext          = 2*(segment_length - 1) + sideSteps + sideSteps;
    int n_test         = n_top + n_bottom + n_right + n_left;
    int n = std::max(n_obs, n_test);
    if (printOutput) printf("n_test:  \t%d\nn_top:  \t%d\nn_right:  \t%d\nn_bottom:\t%d\nn_left:  \t%d\nn_int:  \t%d\nn_ext:  \t%d\nn_obs:  \t%d\nn:      \t%d\n", n_test, n_top, n_right, n_bottom, n_left, n_int, n_ext, n_obs, n);
    this->segments.reserve(num_segments);
    for (int i = 0; i < num_segments; i++) {
        segments.emplace_back(n_obs, n_int, n_ext, n_test, n_top, n_right, n_bottom, n_left, segment_length, deviceComputation);;
    }

}


BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points, bool deviceComputation, double *x_obs, double*y_obs, int n_obs, bool printOutput) {

    // Set variables
    this->protein_structure = protein_structure;
    this->num_segments = num_segments;
    this->total_grid_points = total_grid_points;
    this->deviceComputation = deviceComputation;
    this->printOutput       = printOutput;

    // Setup observation points
    this->n_obs = n_obs;
    this->x_obs = RealMatrix(n_obs);
    this->y_obs = RealMatrix(n_obs);
    for (int i = 0; i < n_obs; i++) this->x_obs.setHostValue(i, x_obs[i]);
    for (int i = 0; i < n_obs; i++) this->y_obs.setHostValue(i, y_obs[i]);
    this->x_obs.toDevice();
    this->y_obs.toDevice();

    // Allocate fields
    E_scat = Field(n_obs);
    H_scat = Field(n_obs);
    E_int = Field(n_obs);
    H_int = Field(n_obs);
    E_inc  = Field(n_obs);
    H_inc  = Field(n_obs);
    for (int i = 0; i < 2; i++) {
        E_scat_pol[i] = Field(n_obs);
        H_scat_pol[i] = Field(n_obs);
        E_int_pol[i]  = Field(n_obs);
        H_int_pol[i]  = Field(n_obs);
        E_inc_pol[i]  = Field(n_obs);
        H_inc_pol[i]  = Field(n_obs);
    }

    // Allocate segments
    int sideSteps      = 0.75*total_grid_points;
    int segment_length = total_grid_points / num_segments;
    int n_top          = segment_length - 2;
    int n_bottom       = segment_length - 2;
    int n_right        = sideSteps - 1;
    int n_left         = sideSteps - 1;
    int n_int          = n_top + n_right + n_bottom + n_left - 16;
    int n_ext          = 2*(segment_length - 1) + sideSteps + sideSteps;
    int n_test         = n_top + n_bottom + n_right + n_left;
    int n = std::max(n_obs, n_test);
    if (printOutput) printf("n_test:  \t%d\nn_top:  \t%d\nn_right:  \t%d\nn_bottom:\t%d\nn_left:  \t%d\nn_int:  \t%d\nn_ext:  \t%d\nn_obs:  \t%d\nn:      \t%d\n", n_test, n_top, n_right, n_bottom, n_left, n_int, n_ext, n_obs, n);
    this->segments.reserve(num_segments);
    for (int i = 0; i < num_segments; i++) {
        segments.emplace_back(n_obs, n_int, n_ext, n_test, n_top, n_right, n_bottom, n_left, segment_length, deviceComputation);;
    }

}

}