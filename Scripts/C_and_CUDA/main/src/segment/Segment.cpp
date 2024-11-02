#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
#include "../../lib/ComplexMatrix.h"
using namespace std;

Segment::Segment() {
    // Empty constructor
}
    
void Segment::free() {
    // Frees all allocated arrays

    x_int.free();           y_int.free();
    x_ext.free();           y_ext.free();
    x_test.free();          y_test.free();
    n_x.free();             n_y.free();
    C.free();
    D.free();
    A.free();
    b.free();
    
}

void Segment::allocate() {
    // Allocates arrays
    
    x_int = RealMatrix(n_int);      y_int = RealMatrix(n_int);
    x_ext = RealMatrix(n_ext);      y_ext = RealMatrix(n_ext);
    x_test = RealMatrix(n_test);    y_test = RealMatrix(n_test);
    n_x = RealMatrix(n_test);       n_y = RealMatrix(n_test);
    C = ComplexMatrix(n_int);          D = ComplexMatrix(n_ext);
    A = RealMatrix(4 * n_test, 2*(n_ext + n_int));
    // Setup right-hand side in linear system
    b = RealMatrix(4 * n_test);

    

}

}