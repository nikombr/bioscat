#ifndef _SEGMENT_H
#define _SEGMENT_H
extern "C" {
#include "RealMatrix.h"
#include "ComplexMatrix.h"
#include "Nanostructure.h"
#include "Field.h"
#include "Constants.h"

class Segment {
    public:
        RealMatrix x_int;
        RealMatrix y_int;
        RealMatrix x_ext;
        RealMatrix y_ext;
        //RealMatrix x_test_top;
        //RealMatrix y_test_top;
        //RealMatrix x_test_right;
        //RealMatrix y_test_right;
        //RealMatrix x_test_bottom;
        //RealMatrix y_test_bottom;
        //RealMatrix x_test_left;
        //RealMatrix y_test_left;
        RealMatrix x_test;
        RealMatrix y_test;
        RealMatrix n_x;
        RealMatrix n_y;
        int num_test_points = 0;
        //int num_interior_points = 0;
        //int num_exterior_points = 0;
        ComplexMatrix C; // Partial solution of the linear system
        ComplexMatrix D; // Partial solution of the linear system
        RealMatrix A; // Linear system matrix in Ax=b
        RealMatrix b; // Linear system vector in Ax=b
        Field E_scat_matrix;
        Field H_scat_matrix;
        Field E_int_matrix; 
        Field H_int_matrix; 
        Field E_inc_vector;
        Field H_inc_vector;
        Field E_ref_vector;
        Field H_ref_vector;
        int polarisation = 1; // Either 1 or 2
        Constants constants;
        //int n_top, n_right, n_bottom, n_left, n_ext, n_int, n_test;
        int n_ext, n_int, n_test;
        int minNumSteps = 10; // Minimum number of steps for sides of segment
        bool deviceComputation;
        int current_segment;



        Segment(); // Empty constructor
        void allocate(); // Allocation of matrices
        void free(); // Free matrices
        void setup(Nanostructure nanostructure, int total_grid_points, int num_segments); // Setup segment
        void computeIncidentFieldVectors(RealMatrix y); // Computes vectors in observation points
        void computeReflectedFieldVectors(RealMatrix y); // Computes vectors in observation points
        void computeScatteredFieldMatrices(RealMatrix x, RealMatrix y, bool far_field_approximation);
        void computeInteriorFieldMatrices(RealMatrix x, RealMatrix y);
        void computeTotalFields();
        void setupRightHandSide();
        void computeFieldsForLinearSystem(); // Computes vectors and matrices in test points
        void setupSystemSubMatrices();
        void solveLinearSystem();
        void setupSystemMatrix();
};
}

#endif