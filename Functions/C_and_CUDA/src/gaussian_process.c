#include <stdlib.h>
#include <stdio.h>
#include "../lib/gaussian_process_inner.h"
#include <cuda_runtime_api.h>

void gaussian_process(double * x, double * y, int n, double * hyper, int num, int dim) {

    gaussian_process_inner(x, y, n, hyper, num, dim);

    printf("Hej fra C\n");

}