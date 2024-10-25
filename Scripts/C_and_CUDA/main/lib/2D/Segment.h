#ifndef _SEGMENT_H
#define _SEGMENT_H
extern "C" {
#include "../RealMatrix.h"
#include "Nanostructure.h"

class Segment {
    public:
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
        RealMatrix n_x;
        RealMatrix n_y;

        Segment(); // Empty constructor
        void allocate(int n_top, int n_right, int n_bottom, int n_left, int n_int, int n_ext); // Allocation of matrices
        void free(); // Free matrices
        void setup(Nanostructure nanostructure, int current_segment, int total_grid_points, int num_segments); // Setup segment
};
}

#endif