#ifndef _SEGMENT_H
#define _SEGMENT_H
#include <cusolverDn.h>
#include "utils/RealMatrix.h"
#include "utils/ComplexMatrix.h"
#include "Nanostructure.h"
#include "Field.h"
#include "utils/Constants.h"
#include "Coordinates.h"
extern "C" {


class Segment {
    public:
        Coordinates aux_int; // Interior auxilliary points
        Coordinates aux_ext; // Exterior auxilliary points
        Coordinates test_points; // Test points
        Coordinates normal_vectors; // Normal vectors of test points
        ComplexMatrix C; // Partial solution of the linear system
        ComplexMatrix D; // Partial solution of the linear system
        RealMatrix A; // Linear system matrix in Ax = b
        RealMatrix b; // Linear system vector in Ax = b
        ComplexMatrix F; // Far field pattern
        Field E_scat_matrix;
        Field H_scat_matrix;
        Field E_int_matrix; 
        Field H_int_matrix; 
        Field E_inc;
        Field H_inc;
        Field E_scat;
        Field H_scat;
        Field E_int; 
        Field H_int;
        int polarisation = 1; // Either 1 or 2
        Constants constants;
        int n_ext, n_int, n_test, n_obs, n_top, n_right, n_bottom, n_left; // Number of points
        int segment_length;
        //int minNumSteps; // Minimum number of steps for sides of segment
        bool deviceComputation = false;
        int current_segment;
        bool printOutput = false;
        cusolverDnHandle_t handle;
        double * A_T_d, *x_d;

        Segment();                                                                        // Constructor
        void allocate();                                                                  // Allocation of matrices
        void free();                                                                      // Free segment and everything allocated
        void setup(Nanostructure nanostructure, int total_grid_points, int num_segments); // Setup segment
        void computeIncidentFields(RealMatrix y);                                         // Computes vectors in given points
        void computeScatteredFields();                                                    // Computes vectors in given points chosen when computing the matrices
        void computeInteriorFields();                                                     // Computes vectors in given points chosen when computing the matrices
        void computeFarFieldPattern(RealMatrix phi);                                      // Computes far field pattern for given angles
        void computeScatteredFieldMatrices(RealMatrix x, RealMatrix y);                   // Computes matrices in given and auxilliary points
        void computeInteriorFieldMatrices(RealMatrix x, RealMatrix y);                    // Computes matrices in given and auxilliary points
        //void computeFieldsForLinearSystem();                                            // Computes vectors and matrices in test points
        void setupRightHandSide();                                                        // Sets up right-hand side b of linear system Ax = b
        void setupSystemMatrix();                                                         // Sets up system matrixs A of linear system Ax = b
        void solveLinearSystem();                                                         // Solves linear system
        //void newWavelength(double lambda);
        //void computeScatteredSubFields();
};

}

#endif