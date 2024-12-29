#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "RealMatrix.h"
#include "Coordinates.h"
extern "C" {
using namespace std;

__global__ void computeInteriorPointsKernel(Coordinates aux_int, Coordinates test_points, double alpha, int n_top, int n_right, int n_bottom, int n_left, int n) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n) {
        double xdiff, ydiff, norm;
        int shift;
        
        if (j < n_top - 4) {
            shift = -2;
        }
        else if (j < n_top + n_right - 8) {
            shift = -6;
        }
        else if (j < n_top + n_right + n_bottom - 12) {
            shift = -10;
        }
        else {
            shift = -14;
        }
        xdiff = test_points.x.getDeviceValue(j - 1 - shift) - test_points.x.getDeviceValue(j + 1 - shift); // Central difference
        ydiff = test_points.y.getDeviceValue(j - 1 - shift) - test_points.y.getDeviceValue(j + 1 - shift); // Central difference
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;
        
        aux_int.x.setDeviceValue(j, test_points.x.getDeviceValue(j - shift) - alpha*ydiff);
        aux_int.y.setDeviceValue(j, test_points.y.getDeviceValue(j - shift) + alpha*xdiff);
    }
}

void computeInteriorPointsCPU(Coordinates aux_int, Coordinates test_points, double alpha, int n_top, int n_right, int n_bottom, int n_left, int n) {
    for (int j = 0; j < n; j++) {
        double xdiff, ydiff, norm;
        int shift;
        
        if (j < n_top - 4) {
            shift = -2;
        }
        else if (j < n_top + n_right - 8) {
            shift = -6;
        }
        else if (j < n_top + n_right + n_bottom - 12) {
            shift = -10;
        }
        else {
            shift = -14;
        }
        xdiff = test_points.x.getHostValue(j - 1 - shift) - test_points.x.getHostValue(j + 1 - shift); // Central difference
        ydiff = test_points.y.getHostValue(j - 1 - shift) - test_points.y.getHostValue(j + 1 - shift); // Central difference
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;
        
        aux_int.x.setHostValue(j, test_points.x.getHostValue(j - shift) - alpha*ydiff);
        aux_int.y.setHostValue(j, test_points.y.getHostValue(j - shift) + alpha*xdiff);
    }
}

void computeInteriorPoints(Coordinates aux_int, Coordinates test_points, double alpha, int n_top, int n_right, int n_bottom, int n_left, int n_int, bool deviceComputation, bool printOutput) {
    int n = n_int;
    if (deviceComputation) { // GPU
        if (printOutput) printf("Computing interior points on the GPU.\n");
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((n + dimBlock.x - 1)/dimBlock.x);
        computeInteriorPointsKernel<<<dimGrid, dimBlock>>>(aux_int, test_points, alpha, n_top, n_right, n_bottom, n_left, n); 
        cudaDeviceSynchronize();
    }
    else { // CPU
        if (printOutput) printf("Computing interior points on the CPU.\n");
        computeInteriorPointsCPU(aux_int, test_points, alpha, n_top, n_right, n_bottom, n_left, n); 

        /*for (int j = 0; j < n; j++) {
            double xdiff, ydiff, norm;
            int shift;
            
            if (j < n_top - 4) {
                shift = -2;
            }
            else if (j < n_top + n_right - 8) {
                shift = -6;
            }
            else if (j < n_top + n_right + n_bottom - 12) {
                shift = -10;
            }
            else {
                shift = -14;
            }
            xdiff = test_points.x.getHostValue(j - 1 - shift) - test_points.x.getHostValue(j + 1 - shift); // Central difference
            ydiff = test_points.y.getHostValue(j - 1 - shift) - test_points.y.getHostValue(j + 1 - shift); // Central difference
            
            norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
            xdiff /= norm;
            ydiff /= norm;
            
            aux_int.x.setHostValue(j, test_points.x.getHostValue(j - shift) - alpha*ydiff);
            aux_int.y.setHostValue(j, test_points.x.getHostValue(j - shift) + alpha*xdiff);

        }*/

    }
}


}