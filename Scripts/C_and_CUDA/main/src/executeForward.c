#include <stdlib.h>
#include <stdio.h>
#include "../lib/GP/gaussian_process_inner.h"
#include "../lib/forward.h"
#include <cuda_runtime_api.h>

void executeForward(double *x, double*y, int n ,char* protein_structure, int num_segments) {

    forward(x,y,n,protein_structure, num_segments);

}