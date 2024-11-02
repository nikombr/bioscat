#ifndef _BIOSCAT_H
#define _BIOSCAT_H
extern "C" {
#include "Nanostructure.h"
#include "Segment.h"
#include "Field.h"
#include "RealMatrix.h"

class BioScat {

    private:
        char * protein_structure;    // Type of protein structure. Either "Retinin2x2" or "demoleus2x2"
        Nanostructure nanostructure; // Information stored about the specific nanostructure 
        int total_grid_points = 1000;       // The number of grid points used along each axis
        Segment * segments;
        int num_segments;
        Field E_scat[2];
        Field H_scat[2];
        Field E_int[2];
        Field H_int[2];
        Field E_inc[2];
        Field H_inc[2];
        Field E_ref[2];
        Field H_ref[2];
        RealMatrix x_obs;
        RealMatrix y_obs;
        bool deviceComputation = false; // True if we should compute on the device
        int polarisation; // Polarisation polarisation, 1 or 2
        

    public:
        BioScat(char* protein_structure, int num_segments);
        ~BioScat();
        void getNanostructure();                                        // Set up nanostructure from protein_structure
        void getSegments();
        void getSegments(Nanostructure nanostructure);
        void forwardSolver(int polarisation);
        void inverseSolver();
        void reset(); // De-allocates everything that only needs to be used in one inverse iteration
        void setupObservationPoints(double *x, double*y, int n);
        void computeScatteredFields();
        void computeReflectedFields();
        void computeIncidentFields();
        void dumpFields();
        void computeReflectance();



};
}
#endif

/*
    private:        
        int n;              // number of points for estimating plane
        double *hyper_h;    // hyperparameters on host
        double *hyper_d;    // hyperparameters on device
        int num;            // number of hyperparameters
        int dim;            // dimension of the problem, either 1 for curve or 2 for plane
        double *x_h;        // x coordinates for estimating plane on host
        double *x_d;        // x coordinates for estimating plane on device
        double *y_h;        // y coordinates for estimating plane on host
        double *y_d;        // y coordinates for estimating plane on device
        double **M_d;       // covariance matrix and later lower triangular matrix from Cholesky factorization on device
        int type_covfunc;
    public:
        double **M_h;       // covariance matrix and later lower triangular matrix from Cholesky factorization on host
        double *M_log;      // M_d[0] on device
        bool device;
        double *p_h;        // random vector and later height of plane in location (x,y) on host
        double *p_d;        // random vector and later height of plane in location (x,y) on device
        GaussianProcess(double* x_h, double* y_h, int n, double* hyper, int num, int dim, int dev, int type_covfunc);  // Constructer, sets default values and allocates
        ~GaussianProcess();                                                                 // Destructer
        void covariance_matrix();                                                           // Computes covariance matrix K
        void cholesky();                                                                    // Does cholesky factorization of K to compute L
        void generate_random_vector();                                                      // Generates random vector p
        void realisation();                                                                 // Computes realisation of the Gaussian process from L
};

#endif*/