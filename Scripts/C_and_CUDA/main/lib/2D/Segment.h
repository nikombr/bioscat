#ifndef _SEGMENTS_H
#define _SEGMENTS_H

#include "../RealMatrix.h"

struct Segment {
    RealMatrix x_int;
    RealMatrix y_int;
    RealMatrix x_ext;
    RealMatrix y_ext;
    RealMatrix x_test_top;
    RealMatrix y_test_top;
    RealMatrix x_test_right;
    RealMatrix y_test_right;
    RealMatrix x_test_bottom;
    RealMatrix y_test_bottom;
    RealMatrix x_test_left;
    RealMatrix y_test_left;

    Segment() {
        std::cout << "Empty segment constructor." << std::endl;
    }
    
    // Constructor with initialization
    Segment(int n1, int n2, int n3, int n4, int n_int, int n_ext)
        : x_int(n_int, 1), y_int(n_int, 1),
          x_ext(n_ext, 1), y_ext(n_ext, 1),
          x_test_top(n1, 1), y_test_top(n1, 1),
          x_test_right(n2, 1), y_test_right(n2, 1),
          x_test_bottom(n3, 1), y_test_bottom(n3, 1),
          x_test_left(n4, 1), y_test_left(n4, 1) {
        std::cout << "Non-empty segment constructor." << std::endl;
    }
};


#endif