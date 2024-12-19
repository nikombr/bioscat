#include "../../lib/profileForward.h"

void executeProfileForward(double *x, double*y, int n ,char* protein_structure, int num_segments, int total_grid_points, double beta, double lambda, int deviceComputation_int) {

    profileForward(x, y, n, protein_structure, num_segments,total_grid_points, beta, lambda, deviceComputation_int);

}