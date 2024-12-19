#include <stdlib.h>
#include <stdio.h>
#include "../../lib/profileGaussianProcess.h"

void executeProfileGaussianProcess(double * x, double * y, int n, double * hyper, int num, int dim, int dev, int type_covfunc) {

    profileGaussianProcess(x, y, n, hyper, num, dim, dev, type_covfunc);

}