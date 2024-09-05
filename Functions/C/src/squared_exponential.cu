
#include <math.h>

__device__ double squared_exponential(double* a, double* b, double tau, double ell) {
    double temp = (a[0]-b[0]) * (a[0]-b[0]) + (a[1]-b[1]) * (a[1]-b[1]);
    return tau*tau*exp(-(temp)/(2*ell*ell));

}