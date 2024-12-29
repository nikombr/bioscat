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
        Coordinates aux_int;                // Interior auxilliary points
        Coordinates aux_ext;                // Exterior auxilliary points
        Coordinates aux_ext_temp;           // Helper array to compute exterior auxilliary points
        Coordinates test_points;            // Test points
        Coordinates normal_vectors;         // Normal vectors of test points
        ComplexMatrix C;                    // Partial solution of the linear system
        ComplexMatrix D;                    // Partial solution of the linear system
        RealMatrix A;                       // Linear system matrix in Ax = b
        RealMatrix b;                       // Linear system vector in Ax = b
        ComplexMatrix far_field_pattern;    // Far field pattern
        Field E_scat_matrix;                // Scattered E-field matrix, E_scat = E_scat_matrix * C 
        Field H_scat_matrix;                // Scattered H-field matrix, H_scat = H_scat_matrix * C 
        Field E_int_matrix;                 // Interior E-field matrix, E_int = E_int_matrix * D
        Field H_int_matrix;                 // Interior H-field matrix, H_int = H_int_matrix * D
        Field E_inc;                        // Incident E-field
        Field H_inc;                        // Incident H-field
        Field E_scat;                       // Scattered E-field
        Field H_scat;                       // Scattered H-field
        Field E_int;                        // Interior E-field
        Field H_int;                        // Interior H-field
        int polarisation = 1;               // Either 1 or 2 respectively corresponding to beta = 0 degrees and beta = 90 degrees
        Constants constants;                // Physical constants
        int n_ext;                          // Number of exterior auxiliary sources
        int n_int;                          // Number of interior auxiliary sources
        int n_test;                         // Number of test points
        int n_top;                          // Number of test points on the top of the nanostructure
        int n_right;                        // Number of test points on the right side of the nanostructure
        int n_bottom;                       // Number of test points on the bottom of the nanostructure
        int n_left;                         // Number of test points on the left side of the nanostructure
        int n_obs;                          // Number of observation points
        int segment_length;                 // Length of the nanostructure on the current segment
        bool deviceComputation = false;     // True if we are computing on the device otherwise false
        int current_segment;                // What number segment is this segment
        bool printOutput = false;           // True if we want to print output to the terminal otherwise false
        cusolverDnHandle_t handle;          // Cusolver handle for when we use cusolver
        double * A_T_d, *x_d;               // Helper arrays for solving linear system on the device

        // Functions
        Segment();                                                                        // Default constructor
        Segment(int n_obs, int n_int, int n_ext, int n_test, int n_top, int n_right, \
                int n_bottom, int n_left, int segment_length, bool deviceComputation);    // Constructor, allocates all arrays
        void free();                                                                      // Free segment and everything allocated
        void setup(Nanostructure nanostructure, int total_grid_points, int num_segments); // Setup segment
        void computeIncidentFields(RealMatrix y);                                         // Computes vectors in given points
        void computeScatteredFields();                                                    // Computes vectors in given points chosen when computing the matrices
        void computeInteriorFields();                                                     // Computes vectors in given points chosen when computing the matrices
        void computeFarFieldPattern(RealMatrix phi);                                      // Computes far field pattern for given angles
        void computeScatteredFieldMatrices(RealMatrix x, RealMatrix y);                   // Computes matrices in given and auxilliary points
        void computeInteriorFieldMatrices(RealMatrix x, RealMatrix y);                    // Computes matrices in given and auxilliary points
        void computeFieldsForLinearSystem();                                              // Computes vectors and matrices in test points
        void setupRightHandSide();                                                        // Sets up right-hand side b of linear system Ax = b
        void setupSystemMatrix();                                                         // Sets up system matrixs A of linear system Ax = b
        void solveLinearSystem();                                                         // Solves linear system
        void newWavelength(double lambda) {
            constants.newWavelength(lambda);
        }
};

}

#endif