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
#include <vector>
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
extern "C" {
using namespace std;

void BioScat::allocateSegments() {
    printf("Please remove this function!\n");
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