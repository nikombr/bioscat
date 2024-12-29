#include "Segment.h"

extern "C" {
using namespace std;

void Segment::free() {

    // Frees all allocated arrays
    aux_int.free();
    aux_ext.free();
    aux_ext_temp.free();
    test_points.free();
    normal_vectors.free();
    C.free();
    D.free();
    A.free();
    b.free();
    far_field_pattern.free();
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