#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

Segment::Segment() {
    // Empty constructor
}
    
void Segment::free() {
    // Frees all allocated arrays

    x_int.free();           y_int.free();
    x_ext.free();           y_ext.free();
    x_test_top.free();      y_test_top.free();
    x_test_right.free();    y_test_right.free();
    x_test_bottom.free();   y_test_bottom.free();
    x_test_left.free();     y_test_left.free();
    n_x.free();             n_y.free();

}

void Segment::allocate(int n_top, int n_right, int n_bottom, int n_left, int n_int, int n_ext) {
    // Allocates arrays
    
    x_int = RealMatrix(n_int);              y_int = RealMatrix(n_int);
    x_ext = RealMatrix(n_ext);              y_ext = RealMatrix(n_ext);
    x_test_top = RealMatrix(n_top);         y_test_top = RealMatrix(n_top);
    n_x = RealMatrix(n_top);                n_y = RealMatrix(n_top);
    x_test_right = RealMatrix(n_right);     y_test_right = RealMatrix(n_right);
    x_test_bottom = RealMatrix(n_bottom);   y_test_bottom = RealMatrix(n_bottom);
    x_test_left = RealMatrix(n_left);       y_test_left = RealMatrix(n_left);

}

}