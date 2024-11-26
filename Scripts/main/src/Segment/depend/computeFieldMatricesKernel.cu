
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include "../../../lib/Segment.h"
#include "../../../lib/utils/RealMatrix.h"
extern "C" {
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

__global__ void computeFieldMatricesKernel(ComplexMatrix F1, ComplexMatrix F2, ComplexMatrix F3, RealMatrix x, RealMatrix x_aux, RealMatrix y, RealMatrix y_aux, double const1, double const2, int rows, int cols, double k) {
    int r = threadIdx.x + blockIdx.x * blockDim.x;
    int c = threadIdx.y + blockIdx.y * blockDim.y;
    if (r < rows && c < cols) {
        double abs_aux, xdiff, ydiff, H_real, H_imag;

        // Get data
        xdiff   = x.getDeviceValue(r) - x_aux.getDeviceValue(c);
        ydiff   = y.getDeviceValue(r) - y_aux.getDeviceValue(c);
        abs_aux = std::sqrt(xdiff * xdiff + ydiff * ydiff);

        // Compute first Hankel functions
        H_real = H02_real(k * abs_aux);
        H_imag = H02_imag(k * abs_aux);

        // Compute the first field
        F1.setDeviceRealValue(r, c, const1 * H_real);
        F1.setDeviceImagValue(r, c, const1 * H_imag);

        // Compute second Hankel functions
        H_real = H12_real(k * abs_aux);
        H_imag = H12_imag(k * abs_aux);

        // Compute the second field
        F2.setDeviceRealValue(r, c,   const2 * 1/abs_aux * H_imag * ydiff);
        F2.setDeviceImagValue(r, c, - const2 * 1/abs_aux * H_real * ydiff);

        // Compute the third field
        F3.setDeviceRealValue(r, c, - const2 * 1/abs_aux * H_imag * xdiff);
        F3.setDeviceImagValue(r, c,   const2 * 1/abs_aux * H_real * xdiff);
    }
}

}