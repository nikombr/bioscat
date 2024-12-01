//#include <stdlib.h>
//#include <stdio.h>
//#include <cuda_runtime_api.h>
//#include <iostream>
//#include <string.h>
//#include "../../lib/cuSolver.h"
#include "../../lib/Segment.h"

extern "C" {
using namespace std;

void Segment::free() {

    // Frees all allocated arrays
    aux_int.free();
    aux_ext.free();
    test_points.free();
    normal_vectors.free();
    C.free();
    D.free();
    A.free();
    b.free();
    F.free();
    E_scat_matrix.free();
    H_scat_matrix.free();
    E_int_matrix.free(); 
    H_int_matrix.free(); 
    E_inc.free();
    H_inc.free();
    E_scat.free();
    H_scat.free();
    E_int.free(); 
    H_int.free();
    cudaFree(A_T_d);
    cudaFree(x_d);
    cusolverDnDestroy(handle);
    
}

}