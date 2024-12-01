#include <stdlib.h>
#include <stdio.h>
#include "../../lib/computeGaussianProcess.h"
#include <cuda_runtime_api.h>

void executeComputeGaussianProcess(double * x, double * y, int n, double * hyper, int num, int dim, int dev, int type_covfunc) {

    computeGaussianProcess(x, y, n, hyper, num, dim, dev, type_covfunc);

}