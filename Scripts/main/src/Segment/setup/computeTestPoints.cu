#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "Nanostructure.h"
#include "RealMatrix.h"
#include "kernels.h"
#include "Coordinates.h"
extern "C" {
using namespace std;


void computeTestPoints(Coordinates test_points, Nanostructure nanostructure, int start, int end, double leftStep, double rightStep, int n_top, int n_right, int n_bottom, int n_left, double left_x_value, double right_x_value, bool deviceComputation, bool printOutput) {

    // Remove end points
    start     += 1;
    end       -= 1;

    if (deviceComputation) { // GPU
        if (printOutput) printf("Computing test points on the GPU.\n");

        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((n_top + dimBlock.x - 1)/dimBlock.x);
        int shift = 0;
        setVectorKernel<<<dimGrid, dimBlock>>>(test_points.x,            shift,         n_top, nanostructure.x, start);
        setVectorKernel<<<dimGrid, dimBlock>>>(test_points.y,            shift,         n_top, nanostructure.f, start);

        dimGrid.x = (n_right + dimBlock.x - 1)/dimBlock.x;
        shift += n_top;
        setConstantKernel<<<dimGrid, dimBlock>>>(test_points.x, shift, n_right, right_x_value);
        setReversedKernel<<<dimGrid, dimBlock>>>(test_points.y, shift, n_right, rightStep);

        dimGrid.x = (n_bottom + dimBlock.x - 1)/dimBlock.x;
        shift += n_right;
        setReversedVectorKernel<<<dimGrid, dimBlock>>>(test_points.x, shift, n_bottom, nanostructure.x, start);
        setConstantKernel<<<dimGrid, dimBlock>>>(test_points.y,       shift, n_bottom, 0.0);

        dimGrid.x = (n_left + dimBlock.x - 1)/dimBlock.x;
        shift += n_bottom;
        setConstantKernel<<<dimGrid, dimBlock>>>(test_points.x, shift, n_left, left_x_value);
        setLinearKernel<<<dimGrid, dimBlock>>>(test_points.y,   shift, n_left, leftStep);

        cudaDeviceSynchronize();

    }
    else { // CPU
        if (printOutput) printf("Computing test points on the CPU.\n");
        int shift = 0;
        setVectorCPU(test_points.x,            shift,         n_top, nanostructure.x, start);
        setVectorCPU(test_points.y,            shift,         n_top, nanostructure.f, start);

        shift += n_top;
        setConstantCPU(test_points.x, shift, n_right, right_x_value);
        setReversedCPU(test_points.y, shift, n_right, rightStep);

     
        shift += n_right;
        setReversedVectorCPU(test_points.x, shift, n_bottom, nanostructure.x, start);
        setConstantCPU(test_points.y,       shift, n_bottom, 0.0);

        shift += n_bottom;
        setConstantCPU(test_points.x, shift, n_left, left_x_value);
        setLinearCPU(test_points.y,   shift, n_left, leftStep);

        /*bool host = !deviceComputation;
        bool device = deviceComputation;
        RealMatrix x_test_top       = RealMatrix(n_top,    host, device);
        RealMatrix y_test_top       = RealMatrix(n_top,    host, device);
        RealMatrix x_test_right     = RealMatrix(n_right,  host, device);
        RealMatrix y_test_right     = RealMatrix(n_right,  host, device);
        RealMatrix x_test_bottom    = RealMatrix(n_bottom, host, device);
        RealMatrix y_test_bottom    = RealMatrix(n_bottom, host, device);
        RealMatrix x_test_left      = RealMatrix(n_left,   host, device);
        RealMatrix y_test_left      = RealMatrix(n_left,   host, device);

        // Compute points along each side
        for (int j = 0; j < n_right; j++) {
            x_test_right.setHostValue(j, right_x_value);
            y_test_right.setHostValue(n_right - j - 1, (j+1)*rightStep);
        }
        
        for (int j = 0; j < n_left; j++) {
            x_test_left.setHostValue(j, left_x_value);
            y_test_left.setHostValue(j, (j+1)*leftStep);
        }

        for (int j = start; j < end; j++) {
            
            x_test_top.setHostValue(j - start, nanostructure.x.getHostValue(j));
            y_test_top.setHostValue(j - start, nanostructure.f.getHostValue(j));
            x_test_bottom.setHostValue(end - j - 1, nanostructure.x.getHostValue(j));
            y_test_bottom.setHostValue(j - start, 0.0);
        }

        // Combine points into combined vector
        int shift = 0;
        for (int j = 0; j < n_top;    j++) {
            test_points.x.setHostValue(j + shift,x_test_top.getHostValue(j));
            test_points.y.setHostValue(j + shift,y_test_top.getHostValue(j));
        }
        shift += n_top;
        for (int j = 0; j < n_right;  j++) {
            test_points.x.setHostValue(j + shift,x_test_right.getHostValue(j));
            test_points.y.setHostValue(j + shift,y_test_right.getHostValue(j));
        }
        shift += n_right;
        for (int j = 0; j < n_bottom; j++) {
            test_points.x.setHostValue(j + shift,x_test_bottom.getHostValue(j));
            test_points.y.setHostValue(j + shift,y_test_bottom.getHostValue(j));
        }
        shift += n_bottom;
        for (int j = 0; j < n_left;   j++) {
            test_points.x.setHostValue(j + shift,x_test_left.getHostValue(j));
            test_points.y.setHostValue(j + shift,y_test_left.getHostValue(j));
        }        
        x_test_top.free();      y_test_top.free();
        x_test_right.free();    y_test_right.free();
        x_test_bottom.free();   y_test_bottom.free();
        x_test_left.free();     y_test_left.free();*/
    }
    
}

}