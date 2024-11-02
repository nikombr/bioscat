#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../../lib/Segment.h"
#include "../../../lib/RealMatrix.h"
using namespace std;


void computeExteriorPointsAndNormalVectors(RealMatrix x_ext, RealMatrix y_ext, RealMatrix n_x, RealMatrix n_y, Nanostructure nanostructure, int start, int end, double alpha, bool deviceComputation) {

    // Allocate array for temporary points to compute exterior auxiliary points
    RealMatrix x_temp = RealMatrix(x_ext.rows + 2);
    RealMatrix y_temp = RealMatrix(y_ext.rows + 2);

    // Set minimum number of steps on the sides of the segment
    int minNumSteps = 10;

    // Determine distance between test points
    double step = nanostructure.x.getHostValue(1) - nanostructure.x.getHostValue(0);

    // Get values at end points
    double startvalue  = nanostructure.f.getHostValue(start);
    double endvalue    = nanostructure.f.getHostValue(end - 1);
    double startxvalue = nanostructure.x.getHostValue(start);
    double endxvalue   = nanostructure.x.getHostValue(end - 1);

    // Determine steps for sides of segment
    int startnum = max(minNumSteps, (int) ceil(startvalue/step));
    int endnum   = max(minNumSteps, (int) ceil(endvalue/step));
    double startstep = startvalue/startnum;
    double endstep   = endvalue/endnum;

    // Get dimensions of each side
    int n_top       = end - start - 2;
    int n_bottom    = end - start - 2;
    int n_right     = endnum - 1;
    int n_left      = startnum - 1;

    if (deviceComputation) { // GPU
        printf("Computing exterior test points on the GPU.\n");

    }
    else { // CPU
        printf("Computing exterior test points on the CPU.\n");

        for (int j = 0; j < endnum; j++) {
            x_temp.setHostValue(end - 1 - start + j + 1, endxvalue);
            y_temp.setHostValue(end - 1 - start + endnum - 1 - j + 1, (j+1)*endstep);
        }
        
        for (int j = 0; j < startnum + 1; j++) {
            x_temp.setHostValue(2*(end - start - 1) + endnum + j + 1, startxvalue);
            y_temp.setHostValue(2*(end - start - 1) + endnum + j + 1, j*startstep);
        }

        x_temp.setHostValue(0, startxvalue);
        y_temp.setHostValue(0, (startnum - 1)*startstep);

        for (int j = start; j < end - 1; j++) {
            
            x_temp.setHostValue(j - start + 1, nanostructure.x.getHostValue(j));
            x_temp.setHostValue(end - start - 1 + endnum + end - j - 2 + 1, nanostructure.x.getHostValue(j+1));
            y_temp.setHostValue(j - start + 1, nanostructure.f.getHostValue(j));
            y_temp.setHostValue(end - start - 1 + endnum + j - start + 1 - 1 + 1, 0.0);
        }

        for (int j = 1; j < x_ext.rows + 1; j++) {
            double xdiff, ydiff, norm;

            xdiff = x_temp.getHostValue(j - 1) - x_temp.getHostValue(j + 1); // Central difference
            ydiff = y_temp.getHostValue(j - 1) - y_temp.getHostValue(j + 1); // Central difference
            
            norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
            xdiff /= norm;
            ydiff /= norm;

            x_ext.setHostValue(j - 1, x_temp.getHostValue(j) + alpha*ydiff);
            y_ext.setHostValue(j - 1, y_temp.getHostValue(j) - alpha*xdiff);

            // Top: Set normal vectors
            if (j >= 2 && j < n_top + 2) {
                n_x.setHostValue(j - 2,   ydiff);
                n_y.setHostValue(j - 2, - xdiff);
            }
        }
        
        // Right side: Set normal vectors
        for (int j = n_top; j < n_top + n_right; j++) {
            
            n_x.setHostValue(j, 1.0);
            n_y.setHostValue(j, 0.0);
        }

        // Bottom: Set normal vectors
        for (int j = n_top + n_right; j < n_top + n_right + n_bottom; j++) {
            n_x.setHostValue(j, 0.0);
            n_y.setHostValue(j, 1.0);
        }

        // Left side: Set normal vectors
        for (int j = n_top + n_right + n_bottom; j < n_top + n_right + n_bottom + n_left; j++) {
            n_x.setHostValue(j, 1.0);
            n_y.setHostValue(j, 0.0);
        }
    }
}


}