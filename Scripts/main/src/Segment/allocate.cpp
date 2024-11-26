//#include <stdlib.h>
//#include <stdio.h>
//#include <cuda_runtime_api.h>
//#include <iostream>
//#include <string.h>
//#include "../../lib/cuSolver.h"
#include "../../lib/Segment.h"

extern "C" {
using namespace std;

void Segment::allocate() {
    // Allocates arrays
    bool host = !deviceComputation;
    bool device = deviceComputation;
    int n = std::max(n_obs, n_test);
    if (printOutput) printf("Allocate: (n_test, n_int, n_ext, n_obs, n) = (%d, %d, %d, %d, %d)\n", n_test, n_int, n_ext, n_obs, n);
    aux_int        = Coordinates(n_int,                       host, device);
    aux_ext        = Coordinates(n_ext,                       host, device);
    test_points    = Coordinates(n_test,                      host, device);
    normal_vectors = Coordinates(n_test,                      host, device);
    C              = ComplexMatrix(n_int,                      host, device);
    D              = ComplexMatrix(n_ext,                      host, device);
    A              = RealMatrix(4 * n_test, 2*(n_ext + n_int), host, device);
    b              = RealMatrix(4 * n_test,                    host, device);
    F              = ComplexMatrix(n,                          host, device);
    E_scat_matrix  = Field(n, n_int,                           host, device);
    H_scat_matrix  = Field(n, n_int,                           host, device);
    E_int_matrix   = Field(n, n_ext,                           host, device);
    H_int_matrix   = Field(n, n_ext,                           host, device);
    E_scat         = Field(n,                                  host, device);
    H_scat         = Field(n,                                  host, device);
    E_int          = Field(n,                                  host, device);
    H_int          = Field(n,                                  host, device);
    E_inc          = Field(n,                                  host, device);
    H_inc          = Field(n,                                  host, device);
    // Prepare for linear system
    cusolverDnCreate(&handle);
    cudaMalloc((void **) &A_T_d,    A.rows * A.cols * sizeof(double));
    cudaMalloc((void **) &x_d,      A.cols * sizeof(double));

}
}