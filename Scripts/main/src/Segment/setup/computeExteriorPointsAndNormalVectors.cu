#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "Nanostructure.h"
#include "Coordinates.h"
#include "RealMatrix.h"
#include "kernels.h"
using namespace std;

__global__ void computeExteriorPointsAndNormalVectorsKernel(Coordinates aux_ext, Coordinates normal_vectors, RealMatrix x_temp, RealMatrix y_temp, int n, double alpha, int n_top) {
    int j = threadIdx.x + blockIdx.x * blockDim.x + 1;
    if (j < n + 1) {
        double xdiff, ydiff, norm;

        xdiff = x_temp.getDeviceValue(j - 1) - x_temp.getDeviceValue(j + 1); // Central difference
        ydiff = y_temp.getDeviceValue(j - 1) - y_temp.getDeviceValue(j + 1); // Central difference
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;

        aux_ext.x.setDeviceValue(j - 1, x_temp.getDeviceValue(j) + alpha*ydiff);
        aux_ext.y.setDeviceValue(j - 1, y_temp.getDeviceValue(j) - alpha*xdiff);

        // Top: Set normal vectors
        if (j >= 2 && j < n_top + 2) {
            normal_vectors.x.setDeviceValue(j - 2,   ydiff);
            normal_vectors.y.setDeviceValue(j - 2, - xdiff);
        }
    }
}

void computeExteriorPointsAndNormalVectorsCPU(Coordinates aux_ext, Coordinates normal_vectors, RealMatrix x_temp, RealMatrix y_temp, int n, double alpha, int n_top) {

    for (int j = 1; j < n + 1; j++) {
        double xdiff, ydiff, norm;

        xdiff = x_temp.getHostValue(j - 1) - x_temp.getHostValue(j + 1); // Central difference
        ydiff = y_temp.getHostValue(j - 1) - y_temp.getHostValue(j + 1); // Central difference
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;

        aux_ext.x.setHostValue(j - 1, x_temp.getHostValue(j) + alpha*ydiff);
        aux_ext.y.setHostValue(j - 1, y_temp.getHostValue(j) - alpha*xdiff);

        // Top: Set normal vectors
        if (j >= 2 && j < n_top + 2) {
            normal_vectors.x.setHostValue(j - 2,   ydiff);
            normal_vectors.y.setHostValue(j - 2, - xdiff);
        }
    }
}


void computeExteriorPointsAndNormalVectors(Coordinates aux_ext, Coordinates normal_vectors, Coordinates aux_ext_temp, Nanostructure nanostructure, int start, int end, double alpha, double leftStep, double rightStep, int leftNum, int rightNum, int n_top, int n_right, int n_bottom, int n_left, double left_x_value, double right_x_value, int n_ext, bool deviceComputation, bool printOutput) {
   

    if (deviceComputation) { // GPU
        if (printOutput) printf("Computing exterior test points on the GPU.\n");

        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((rightNum + dimBlock.x - 1)/dimBlock.x);

        setConstantKernel<<<dimGrid, dimBlock>>>(aux_ext_temp.x, end - start, rightNum, right_x_value);
        setReversedKernel<<<dimGrid, dimBlock>>>(aux_ext_temp.y, end - start, rightNum, rightStep);

        dimGrid.x = (leftNum + dimBlock.x - 1)/dimBlock.x;
        setConstantKernel<<<dimGrid, dimBlock>>>(aux_ext_temp.x, 2*(end - start - 1) + rightNum + 1, leftNum, left_x_value);
        setLinearKernel2<<<dimGrid, dimBlock>>>(aux_ext_temp.y,  2*(end - start - 1) + rightNum + 1, leftNum, leftStep);

        setConstantKernel<<<1, 1>>>(aux_ext_temp.x, 0, 1, left_x_value);
        setConstantKernel<<<1, 1>>>(aux_ext_temp.y, 0, 1, (leftNum - 1)*leftStep);

        int n = end - 1 - start;
        dimGrid.x = (n + dimBlock.x - 1)/dimBlock.x;
        setVectorKernel<<<dimGrid, dimBlock>>>(        aux_ext_temp.x, 1,                    n, nanostructure.x, start);
        setReversedVectorKernel<<<dimGrid, dimBlock>>>(aux_ext_temp.x, end - start + rightNum, n, nanostructure.x, start+1);
        setVectorKernel<<<dimGrid, dimBlock>>>(        aux_ext_temp.y, 1,                    n, nanostructure.f, start);
        setConstantKernel<<<dimGrid, dimBlock>>>(      aux_ext_temp.y, end - start + rightNum, n, 0.0);

        cudaDeviceSynchronize();

        dimGrid.x = (n_ext + dimBlock.x - 1)/dimBlock.x;
        computeExteriorPointsAndNormalVectorsKernel<<<dimGrid, dimBlock>>>(aux_ext, normal_vectors, aux_ext_temp.x, aux_ext_temp.y, n_ext, alpha, n_top);

        // Set normal vectors
        dimGrid.x = (n_right + dimBlock.x - 1)/dimBlock.x;
        setConstantKernel<<<dimGrid, dimBlock>>>(normal_vectors.x, n_top, n_right, 1.0);
        setConstantKernel<<<dimGrid, dimBlock>>>(normal_vectors.y, n_top, n_right, 0.0);

        dimGrid.x = (n_bottom + dimBlock.x - 1)/dimBlock.x;
        setConstantKernel<<<dimGrid, dimBlock>>>(normal_vectors.x, n_top + n_right, n_bottom, 0.0);
        setConstantKernel<<<dimGrid, dimBlock>>>(normal_vectors.y, n_top + n_right, n_bottom, 1.0);

        dimGrid.x = (n_left + dimBlock.x - 1)/dimBlock.x;
        setConstantKernel<<<dimGrid, dimBlock>>>(normal_vectors.x, n_top + n_right + n_bottom, n_left, 1.0);
        setConstantKernel<<<dimGrid, dimBlock>>>(normal_vectors.y, n_top + n_right + n_bottom, n_left, 0.0);

        cudaDeviceSynchronize();
    }
    else { // CPU
        if (printOutput) printf("Computing exterior test points on the CPU.\n");

        setConstantCPU(aux_ext_temp.x, end - start, rightNum, right_x_value);
        setReversedCPU(aux_ext_temp.y, end - start, rightNum, rightStep);

      
        setConstantCPU(aux_ext_temp.x, 2*(end - start - 1) + rightNum + 1, leftNum, left_x_value);
        setLinearCPU2(aux_ext_temp.y,  2*(end - start - 1) + rightNum + 1, leftNum, leftStep);

        setConstantCPU(aux_ext_temp.x, 0, 1, left_x_value);
        setConstantCPU(aux_ext_temp.y, 0, 1, (leftNum - 1)*leftStep);

        int n = end - 1 - start;
        setVectorCPU(        aux_ext_temp.x, 1,                    n, nanostructure.x, start);
        setReversedVectorCPU(aux_ext_temp.x, end - start + rightNum, n, nanostructure.x, start+1);
        setVectorCPU(        aux_ext_temp.y, 1,                    n, nanostructure.f, start);
        setConstantCPU(      aux_ext_temp.y, end - start + rightNum, n, 0.0);

        computeExteriorPointsAndNormalVectorsCPU(aux_ext, normal_vectors, aux_ext_temp.x, aux_ext_temp.y, n_ext, alpha, n_top);

        // Set normal vectors

        setConstantCPU(normal_vectors.x, n_top, n_right, 1.0);
        setConstantCPU(normal_vectors.y, n_top, n_right, 0.0);

 
        setConstantCPU(normal_vectors.x, n_top + n_right, n_bottom, 0.0);
        setConstantCPU(normal_vectors.y, n_top + n_right, n_bottom, 1.0);

    
        setConstantCPU(normal_vectors.x, n_top + n_right + n_bottom, n_left, 1.0);
        setConstantCPU(normal_vectors.y, n_top + n_right + n_bottom, n_left, 0.0);

        

        /*for (int j = 0; j < rightNum; j++) {
            aux_ext_temp.x.setHostValue(end - start + j, right_x_value);
            aux_ext_temp.y.setHostValue(end - start + rightNum - j - 1, (j+1)*rightStep);
        }
        
        for (int j = 0; j < leftNum + 1; j++) {
            aux_ext_temp.x.setHostValue(2*(end - start - 1) + rightNum + j + 1, left_x_value);
            aux_ext_temp.y.setHostValue(2*(end - start - 1) + rightNum + j + 1, j*leftStep);
        }

        aux_ext_temp.x.setHostValue(0, left_x_value);
        aux_ext_temp.y.setHostValue(0, (leftNum - 1)*leftStep);

        for (int j = start; j < end - 1; j++) {
            
            aux_ext_temp.x.setHostValue(j - start + 1, nanostructure.x.getHostValue(j));
            aux_ext_temp.x.setHostValue(end - start - 1 + rightNum + end - j - 1, nanostructure.x.getHostValue(j+1));
            aux_ext_temp.y.setHostValue(j - start + 1, nanostructure.f.getHostValue(j));
            aux_ext_temp.y.setHostValue(end - start + rightNum + j - start, 0.0);
        }

        for (int j = 1; j < n_ext + 1; j++) {
            double xdiff, ydiff, norm;

            xdiff = aux_ext_temp.x.getHostValue(j - 1) - aux_ext_temp.x.getHostValue(j + 1); // Central difference
            ydiff = aux_ext_temp.y.getHostValue(j - 1) - aux_ext_temp.y.getHostValue(j + 1); // Central difference
            
            norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
            xdiff /= norm;
            ydiff /= norm;

            aux_ext.x.setHostValue(j - 1, aux_ext_temp.x.getHostValue(j) + alpha*ydiff);
            aux_ext.y.setHostValue(j - 1, aux_ext_temp.y.getHostValue(j) - alpha*xdiff);

            // Top: Set normal vectors
            if (j >= 2 && j < n_top + 2) {
                normal_vectors.x.setHostValue(j - 2,   ydiff);
                normal_vectors.y.setHostValue(j - 2, - xdiff);
            }
        }
        
        // Right side: Set normal vectors
        for (int j = n_top; j < n_top + n_right; j++) {
            
            normal_vectors.x.setHostValue(j, 1.0);
            normal_vectors.y.setHostValue(j, 0.0);
        }

        // Bottom: Set normal vectors
        for (int j = n_top + n_right; j < n_top + n_right + n_bottom; j++) {
            normal_vectors.x.setHostValue(j, 0.0);
            normal_vectors.y.setHostValue(j, 1.0);
        }

        // Left side: Set normal vectors
        for (int j = n_top + n_right + n_bottom; j < n_top + n_right + n_bottom + n_left; j++) {
            normal_vectors.x.setHostValue(j, 1.0);
            normal_vectors.y.setHostValue(j, 0.0);
        }*/
    }

    //x_temp.free();
    //y_temp.free();
}


}