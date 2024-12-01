#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../../../lib/Segment.h"
#include "../../../lib/BioScat.h"
#include "../../../lib/utils/RealMatrix.h"
extern "C" {
using namespace std;

__global__ void combinePolarisationKernel(ComplexMatrix combined, ComplexMatrix P1, ComplexMatrix P2, int rows, double cosBeta, double sinBeta) {
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < rows) {
        combined.setDeviceRealValue(i, P1.getDeviceRealValue(i)*cosBeta + P2.getDeviceRealValue(i)*sinBeta);
        combined.setDeviceImagValue(i, P1.getDeviceImagValue(i)*cosBeta + P2.getDeviceImagValue(i)*sinBeta);
    }
}

void combinePolarisation(Field * pol, Field combined, double beta, bool deviceComputation) {

    int rows = combined.x.rows;
    double cosBeta = cos(beta);
    double sinBeta = sin(beta);
    
    if (deviceComputation) {

        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rows + dimBlock.x - 1)/dimBlock.x);

        combinePolarisationKernel<<<dimGrid, dimBlock>>>(combined.x, pol[0].x, pol[1].x, rows, cosBeta, sinBeta);
        combinePolarisationKernel<<<dimGrid, dimBlock>>>(combined.y, pol[0].y, pol[1].y, rows, cosBeta, sinBeta);
        combinePolarisationKernel<<<dimGrid, dimBlock>>>(combined.z, pol[0].z, pol[1].z, rows, cosBeta, sinBeta);
        
        cudaDeviceSynchronize();
        //combined.toHost();

    }
    else {
        #pragma omp parallel for
        for (int i = 0; i < rows; i++) {           
            combined.x.setHostRealValue(i, pol[0].x.getHostRealValue(i)*cosBeta + pol[1].x.getHostRealValue(i)*sinBeta);
            combined.y.setHostRealValue(i, pol[0].y.getHostRealValue(i)*cosBeta + pol[1].y.getHostRealValue(i)*sinBeta);
            combined.z.setHostRealValue(i, pol[0].z.getHostRealValue(i)*cosBeta + pol[1].z.getHostRealValue(i)*sinBeta);
            combined.x.setHostImagValue(i, pol[0].x.getHostImagValue(i)*cosBeta + pol[1].x.getHostImagValue(i)*sinBeta);
            combined.y.setHostImagValue(i, pol[0].y.getHostImagValue(i)*cosBeta + pol[1].y.getHostImagValue(i)*sinBeta);
            combined.z.setHostImagValue(i, pol[0].z.getHostImagValue(i)*cosBeta + pol[1].z.getHostImagValue(i)*sinBeta);  
        }
    }

}

}