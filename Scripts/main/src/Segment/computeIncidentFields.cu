#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "../../lib/Segment.h"
extern "C" {
#include "../../lib/utils/RealMatrix.h"
using namespace std;

__global__ void computeFieldKernelReal(ComplexMatrix field, RealMatrix y, int n, double interior_constant, double exterior_constant) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n) {
        field.setDeviceRealValue(j, exterior_constant * cos(interior_constant * y.getDeviceValue(j)));
    }
}

__global__ void computeFieldKernelImag(ComplexMatrix field, RealMatrix y, int n, double interior_constant, double exterior_constant) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n) {
        field.setDeviceImagValue(j, exterior_constant * sin(interior_constant * y.getDeviceValue(j)));
    }
}


void Segment::computeIncidentFields(RealMatrix y) {
    
    int rows = y.rows;

    if (deviceComputation) {
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x);
        if (polarisation == 1) {
            computeFieldKernelReal<<<dimGrid, dimBlock>>>(E_inc.z, y, rows, constants.k0, 1.0);
            computeFieldKernelImag<<<dimGrid, dimBlock>>>(E_inc.z, y, rows, constants.k0, 1.0);
            computeFieldKernelReal<<<dimGrid, dimBlock>>>(H_inc.x, y, rows, constants.k0, -1/constants.eta0);
            computeFieldKernelImag<<<dimGrid, dimBlock>>>(H_inc.x, y, rows, constants.k0, -1/constants.eta0);
        }
        else if (polarisation == 2) {
            computeFieldKernelReal<<<dimGrid, dimBlock>>>(E_inc.x, y, rows, constants.k0, 1.0);
            computeFieldKernelImag<<<dimGrid, dimBlock>>>(E_inc.x, y, rows, constants.k0, 1.0);
            computeFieldKernelReal<<<dimGrid, dimBlock>>>(H_inc.z, y, rows, constants.k0, 1/constants.eta0);
            computeFieldKernelImag<<<dimGrid, dimBlock>>>(H_inc.z, y, rows, constants.k0, 1/constants.eta0);
        }

        cudaDeviceSynchronize();
    }
    else {
        if (polarisation == 1) {
            for (int j = 0; j < rows; j++) E_inc.z.setHostRealValue(j,                     cos(constants.k0 * y.getHostValue(j)));
            for (int j = 0; j < rows; j++) E_inc.z.setHostImagValue(j,                     sin(constants.k0 * y.getHostValue(j)));
            for (int j = 0; j < rows; j++) H_inc.x.setHostRealValue(j, -1/constants.eta0 * cos(constants.k0 * y.getHostValue(j)));
            for (int j = 0; j < rows; j++) H_inc.x.setHostImagValue(j, -1/constants.eta0 * sin(constants.k0 * y.getHostValue(j)));
        }
        else {
            for (int j = 0; j < rows; j++) E_inc.x.setHostRealValue(j,                    cos(constants.k0 * y.getHostValue(j)));
            for (int j = 0; j < rows; j++) E_inc.x.setHostImagValue(j,                    sin(constants.k0 * y.getHostValue(j)));
            for (int j = 0; j < rows; j++) H_inc.z.setHostRealValue(j, 1/constants.eta0 * cos(constants.k0 * y.getHostValue(j)));
            for (int j = 0; j < rows; j++) H_inc.z.setHostImagValue(j, 1/constants.eta0 * sin(constants.k0 * y.getHostValue(j)));
        }
    }
}


}