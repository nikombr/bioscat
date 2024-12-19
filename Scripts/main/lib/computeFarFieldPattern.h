#ifndef _COMPUTE_FAR_FIELD_PATTERN_H
#define _COMPUTE_FAR_FIELD_PATTERN_H

void computeFarFieldPattern(double * phi, int n ,char * protein_structure, int num_segments, int total_grid_points, double beta, double lambda, int deviceComputation_int, int printOutput_int);


#endif