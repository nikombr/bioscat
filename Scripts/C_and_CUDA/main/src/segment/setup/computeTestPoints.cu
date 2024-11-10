#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../../lib/Segment.h"
#include "../../../lib/RealMatrix.h"
using namespace std;


void computeTestPoints(RealMatrix x_test, RealMatrix y_test, Nanostructure nanostructure, int start, int end, bool deviceComputation) {
    bool printOutput = false;
    // Set minimum number of steps on the sides of the segment
    int minNumSteps = 15;

    // Determine distance between test points
    double step = nanostructure.x.getHostValue(1) - nanostructure.x.getHostValue(0);

    // Get values at end points
    double startvalue  = nanostructure.f.getHostValue(start);
    double endvalue    = nanostructure.f.getHostValue(end - 1);
    double startxvalue = nanostructure.x.getHostValue(start);
    double endxvalue   = nanostructure.x.getHostValue(end - 1);

    // Determine steps for sides of segment
    int startnum = max(minNumSteps, (int) ceil(startvalue/step));
    startnum = minNumSteps;
    int endnum   = max(minNumSteps, (int) ceil(endvalue/step));
    endnum = minNumSteps;
    double startstep = startvalue/startnum;
    double endstep   = endvalue/endnum;

    // Allocate array for test points along each side
    int n_top       = end - start - 2;
    int n_bottom    = end - start - 2;
    int n_right     = endnum - 1;
    int n_left      = startnum - 1;
    RealMatrix x_test_top = RealMatrix(n_top);
    RealMatrix y_test_top = RealMatrix(n_top);
    RealMatrix x_test_right = RealMatrix(n_right);
    RealMatrix y_test_right = RealMatrix(n_right);
    RealMatrix x_test_bottom = RealMatrix(n_bottom);
    RealMatrix y_test_bottom = RealMatrix(n_bottom);
    RealMatrix x_test_left = RealMatrix(n_left);
    RealMatrix y_test_left = RealMatrix(n_left);

    // Remove end points
    start     += 1;
    end       -= 1;

    if (deviceComputation) { // GPU
        if (printOutput) printf("Computing test points on the GPU.\n");

    }
    else { // CPU
        if (printOutput) printf("Computing test points on the CPU.\n");

        // Compute points along each side
        for (int j = 0; j < n_right; j++) {
            x_test_right.setHostValue(j, endxvalue);
            y_test_right.setHostValue(n_right - j - 1, (j+1)*endstep);
        }
        
        for (int j = 0; j < n_left; j++) {
            x_test_left.setHostValue(j, startxvalue);
            y_test_left.setHostValue(j, (j+1)*startstep);
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
            x_test.setHostValue(j + shift,x_test_top.getHostValue(j));
            y_test.setHostValue(j + shift,y_test_top.getHostValue(j));
        }
        shift += n_top;
        for (int j = 0; j < n_right;  j++) {
            x_test.setHostValue(j + shift,x_test_right.getHostValue(j));
            y_test.setHostValue(j + shift,y_test_right.getHostValue(j));
        }
        shift += n_right;
        for (int j = 0; j < n_bottom; j++) {
            x_test.setHostValue(j + shift,x_test_bottom.getHostValue(j));
            y_test.setHostValue(j + shift,y_test_bottom.getHostValue(j));
        }
        shift += n_bottom;
        for (int j = 0; j < n_left;   j++) {
            x_test.setHostValue(j + shift,x_test_left.getHostValue(j));
            y_test.setHostValue(j + shift,y_test_left.getHostValue(j));
        }        

    }
    x_test_top.free();      y_test_top.free();
    x_test_right.free();    y_test_right.free();
    x_test_bottom.free();   y_test_bottom.free();
    x_test_left.free();     y_test_left.free();
}

}