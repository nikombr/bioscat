#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

__host__ __device__ double H02_real(double x) {
    // Computes real part of Hankel function of order zero and second kind
    int n = 0;
    double Jn = jn(n, x); // Compute Bessel functions of the first (Jn) 
    return Jn;
}

__host__ __device__ double H02_imag(double x) {
    // Computes imaginary part of Hankel function of order zero and second kind
    int n = 0;
    double Yn = yn(n, x); // Compute Bessel functions of the second (Yn) kind
    return -Yn;
}


__host__ __device__ double H12_real(double x) {
    // Computes real part of Hankel function of order one and second kind
    int n = 1;
    double Jn = jn(n, x); // Compute Bessel functions of the first (Jn) 
    return Jn; 
}


 __host__ __device__ double H12_imag(double x) {
    // Computes imaginary part of Hankel function of order one and second kind
    int n = 1;
    double Yn = yn(n, x); // Compute Bessel functions of the second (Yn) kind
    return -Yn;
}

void Segment::computeScatteredFieldMatrices(RealMatrix x, RealMatrix y, bool far_field_approximation) {
    
    int rows = y.rows;
    int cols = y_int.rows;
    

    // Determine which matrices to be allocated
    bool Ex_bool = polarisation == 1 ? false : true;
    bool Ey_bool = polarisation == 1 ? false : true;
    bool Ez_bool = polarisation == 1 ? true  : false;
    bool Hx_bool = polarisation == 1 ? true  : false;
    bool Hy_bool = polarisation == 1 ? true  : false;
    bool Hz_bool = polarisation == 1 ? false : true;

    E_scat_matrix = Field(rows, cols, Ex_bool, Ey_bool, Ez_bool);
    H_scat_matrix = Field(rows, cols, Hx_bool, Hy_bool, Hz_bool);

    if (polarisation == 1) {
        #pragma omp parallel for collapse(2) 
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double abs_int, abs_int_ref, xdiff, ydiff, ydiff_ref, H_real, H_imag, H_real_ref, H_imag_ref, val;

                // Get data
                xdiff       = x.getHostValue(r) - x_int.getHostValue(c);
                ydiff       = y.getHostValue(r) - y_int.getHostValue(c);
                ydiff_ref   = y.getHostValue(r) + y_int.getHostValue(c);
                abs_int     = std::sqrt(xdiff*xdiff + ydiff*ydiff);
                abs_int_ref = std::sqrt(xdiff*xdiff + ydiff_ref*ydiff_ref);

                // Compute first Hankel functions
                H_real     = H02_real(constants.k0*abs_int);
                H_real_ref = H02_real(constants.k0*abs_int_ref);
                H_imag     = H02_imag(constants.k0*abs_int);
                H_imag_ref = H02_imag(constants.k0*abs_int_ref);
                
                val = H_real + constants.Gamma_ref * H_real_ref;
                E_scat_matrix.z.setHostRealValue(r, c, val);
                val = H_imag + constants.Gamma_ref * H_imag_ref;
                E_scat_matrix.z.setHostImagValue(r, c, val);

                // Compute second Hankel functions
                H_real     = H12_real(constants.k0*abs_int);
                H_real_ref = H12_real(constants.k0*abs_int_ref);
                H_imag     = H12_imag(constants.k0*abs_int);
                H_imag_ref = H12_imag(constants.k0*abs_int_ref);

                val = 1/constants.eta0 * (1/abs_int      * H_imag     * ydiff + \
                     constants.Gamma_ref * 1/abs_int_ref * H_imag_ref * ydiff_ref);
                H_scat_matrix.x.setHostRealValue(r, c, val);
                val = -1/constants.eta0 * (1/abs_int     * H_real     * ydiff + \
                     constants.Gamma_ref * 1/abs_int_ref * H_real_ref * ydiff_ref);
                H_scat_matrix.x.setHostImagValue(r, c, val);

                val = -1/constants.eta0 * xdiff * (1/abs_int      * H_imag      + \
                             constants.Gamma_ref * 1/abs_int_ref  * H_imag_ref);
                H_scat_matrix.y.setHostRealValue(r, c, val);
                val = 1/constants.eta0 * xdiff * (1/abs_int     * H_real      + \
                            constants.Gamma_ref * 1/abs_int_ref * H_real_ref);
                H_scat_matrix.y.setHostImagValue(r, c, val);
            }
        }

    }
    else if (polarisation == 2) {

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double abs_int, abs_int_ref, xdiff, ydiff, ydiff_ref, H_real, H_imag, H_real_ref, H_imag_ref, val;

                // Get data
                xdiff       = x.getHostValue(r) - x_int.getHostValue(c);
                ydiff       = y.getHostValue(r) - y_int.getHostValue(c);
                ydiff_ref   = y.getHostValue(r) + y_int.getHostValue(c);
                abs_int     = std::sqrt(xdiff*xdiff + ydiff*ydiff);
                abs_int_ref = std::sqrt(xdiff*xdiff + ydiff_ref*ydiff_ref);

                // Compute first Hankel functions
                H_real     = H02_real(constants.k0*abs_int);
                H_real_ref = H02_real(constants.k0*abs_int_ref);
                H_imag     = H02_imag(constants.k0*abs_int);
                H_imag_ref = H02_imag(constants.k0*abs_int_ref);
                
                val = H_real + constants.Gamma_ref * H_real_ref;
                H_scat_matrix.z.setHostRealValue(r, c, val);
                val = H_imag + constants.Gamma_ref * H_imag_ref;
                H_scat_matrix.z.setHostImagValue(r, c, val);

                // Compute second Hankel functions
                H_real     = H12_real(constants.k0*abs_int);
                H_real_ref = H12_real(constants.k0*abs_int_ref);
                H_imag     = H12_imag(constants.k0*abs_int);
                H_imag_ref = H12_imag(constants.k0*abs_int_ref);

                double constant = constants.mu0/(constants.eta0*constants.epsilon0);

                val = -constant * (1/abs_int      * H_imag     * ydiff + \
                     constants.Gamma_ref * 1/abs_int_ref * H_imag_ref * ydiff_ref);
                E_scat_matrix.x.setHostRealValue(r, c, val);
                val = constant * (1/abs_int     * H_real     * ydiff + \
                     constants.Gamma_ref * 1/abs_int_ref * H_real_ref * ydiff_ref);
                E_scat_matrix.x.setHostImagValue(r, c, val);

                val = constant * xdiff * (1/abs_int      * H_imag      + \
                             constants.Gamma_ref * 1/abs_int_ref  * H_imag_ref);
                E_scat_matrix.y.setHostRealValue(r, c, val);
                val = -constant * xdiff * (1/abs_int     * H_real      + \
                            constants.Gamma_ref * 1/abs_int_ref * H_real_ref);
                E_scat_matrix.y.setHostImagValue(r, c, val);
            }
        }

    }
    else {
        printf("Please input 1 or 2 for the polarisation!\n");
    }

}

void Segment::computeInteriorFieldMatrices(RealMatrix x, RealMatrix y) {
    
    int rows = y.rows;
    int cols = y_ext.rows;
    

    // Determine which matrices to be allocated
    bool Ex_bool = polarisation == 1 ? false : true;
    bool Ey_bool = polarisation == 1 ? false : true;
    bool Ez_bool = polarisation == 1 ? true  : false;
    bool Hx_bool = polarisation == 1 ? true  : false;
    bool Hy_bool = polarisation == 1 ? true  : false;
    bool Hz_bool = polarisation == 1 ? false : true;

    E_int_matrix = Field(rows, cols, Ex_bool, Ey_bool, Ez_bool);
    H_int_matrix = Field(rows, cols, Hx_bool, Hy_bool, Hz_bool);

    if (polarisation == 1) {
        #pragma omp parallel for collapse(2) 
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double abs_ext, xdiff, ydiff, H_real, H_imag, val;

                // Get data
                xdiff   = x.getHostValue(r) - x_ext.getHostValue(c);
                ydiff   = y.getHostValue(r) - y_ext.getHostValue(c);
                abs_ext = std::sqrt(xdiff*xdiff + ydiff*ydiff);

                // Compute first Hankel functions
                H_real = H02_real(constants.k1*abs_ext);
                H_imag = H02_imag(constants.k1*abs_ext);
                
                val = H_real;
                E_int_matrix.z.setHostRealValue(r, c, val);
                val = H_imag;
                E_int_matrix.z.setHostImagValue(r, c, val);

                // Compute second Hankel functions
                H_real = H12_real(constants.k1*abs_ext);
                H_imag = H12_imag(constants.k1*abs_ext);

                val =   constants.n1/constants.eta0 * 1/abs_ext * ydiff * H_imag;
                H_int_matrix.x.setHostRealValue(r, c, val);
                val = - constants.n1/constants.eta0 * 1/abs_ext * ydiff * H_real;
                H_int_matrix.x.setHostImagValue(r, c, val);

                val = -constants.n1/constants.eta0 * 1/abs_ext * xdiff * H_imag;
                H_int_matrix.y.setHostRealValue(r, c, val);
                val =  constants.n1/constants.eta0 * 1/abs_ext * xdiff * H_real;
                H_int_matrix.y.setHostImagValue(r, c, val);
            }
        }
    }
    else if (polarisation == 2) {

        #pragma omp parallel for collapse(2) 
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                double abs_ext, xdiff, ydiff, H_real, H_imag, val;

                // Get data
                xdiff   = x.getHostValue(r) - x_ext.getHostValue(c);
                ydiff   = y.getHostValue(r) - y_ext.getHostValue(c);
                abs_ext = std::sqrt(xdiff*xdiff + ydiff*ydiff);

                // Compute first Hankel functions
                H_real = H02_real(constants.k1*abs_ext);
                H_imag = H02_imag(constants.k1*abs_ext);
                
                val = H_real;
                H_int_matrix.z.setHostRealValue(r, c, val);
                val = H_imag;
                H_int_matrix.z.setHostImagValue(r, c, val);

                // Compute second Hankel functions
                H_real = H12_real(constants.k1*abs_ext);
                H_imag = H12_imag(constants.k1*abs_ext);

                double constant = constants.n1*constants.mu0/(constants.eta0*constants.epsilon0);

                val =  - constant * 1/abs_ext * ydiff * H_imag;
                E_int_matrix.x.setHostRealValue(r, c, val);
                val =  constant * 1/abs_ext * ydiff * H_real;
                E_int_matrix.x.setHostImagValue(r, c, val);

                val = constant * 1/abs_ext * xdiff * H_imag;
                E_int_matrix.y.setHostRealValue(r, c, val);
                val = - constant * 1/abs_ext * xdiff * H_real;
                E_int_matrix.y.setHostImagValue(r, c, val);
            }
        }

    }
    else {
        printf("Please input 1 or 2 for the polarisation!\n");
    }

}


}