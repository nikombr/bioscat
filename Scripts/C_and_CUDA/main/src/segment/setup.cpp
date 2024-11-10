#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <math.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

void computeExteriorPointsAndNormalVectors(RealMatrix x_ext, RealMatrix y_ext, RealMatrix n_x, RealMatrix n_y, Nanostructure nanostructure, int start, int end, double alpha, bool deviceComputation);
void computeTestPoints(RealMatrix x_test, RealMatrix y_test, Nanostructure nanostructure, int start, int end, bool deviceComputation);
void computeInteriorPoints(RealMatrix x_int, RealMatrix y_int, RealMatrix x_test, RealMatrix y_test,double alpha, bool deviceComputation, int n_top, int n_right, int n_bottom, int n_left);

void Segment::setup(Nanostructure nanostructure, int total_grid_points, int num_segments) {
    // Sets up test points and auxiliary points for the segment given a nanostructure
    // that is both on the host and on the device

    // Initialize variables
    int segment_length = total_grid_points / num_segments;
    int start, end, startnum, endnum;
    double startvalue, endvalue, step, startstep, endstep, startxvalue, endxvalue, alpha;

    // Determine distance of auxilliary points from nearest test point
    step = nanostructure.x.getHostValue(1) - nanostructure.x.getHostValue(0);
    //alpha = std::min(10e-8,2*std::min(startstep,endstep));

    // Determine wheree the segment starts and ends
    start = current_segment * segment_length;
    end = min(start + segment_length, total_grid_points);
    //printf("(start, end) = (%d, %d)",start,end);

    // Get values at end points
    startvalue  = nanostructure.f.getHostValue(start);
    endvalue    = nanostructure.f.getHostValue(end - 1);
    startxvalue = nanostructure.x.getHostValue(start);
    endxvalue   = nanostructure.x.getHostValue(end - 1);

    // Determine steps for sides of segment
    startnum = max(minNumSteps, (int) ceil(startvalue/step));
    endnum   = max(minNumSteps, (int) ceil(endvalue/step));
    startnum = minNumSteps;
    endnum = minNumSteps;
    startstep = startvalue/startnum;
    endstep   = endvalue/endnum;
    alpha = std::min((double)2*step,(double)1.0e-9),2*std::min(startstep,endstep);
    
    // Allocate arrays
    int n_top       = end - start - 2;
    int n_bottom    = end - start - 2;
    int n_right     = endnum - 1;
    int n_left      = startnum - 1;
    n_int           = n_top + n_right + n_bottom + n_left - 16;
    n_ext           = 2*(end - start - 1) + endnum + startnum;
    n_test          = n_top + n_bottom + n_right + n_left;
    allocate();

    // Compute exterior points
    computeExteriorPointsAndNormalVectors(x_ext, y_ext, n_x, n_y, nanostructure, start, end, alpha, deviceComputation);

    // Compute test points and temporary test points for interior points
    computeTestPoints(x_test, y_test, nanostructure, start, end, deviceComputation);
 
    // Compute interior points
    computeInteriorPoints(x_int, y_int, x_test, y_test, alpha, deviceComputation, n_top, n_right,  n_bottom,  n_left);

    x_ext.toDevice();
    y_ext.toDevice();
    x_int.toDevice();
    y_int.toDevice();
    x_test.toDevice();
    y_test.toDevice();
    n_x.toDevice();
    n_y.toDevice();
    /*
    for (int j = 0; j < n_int; j++) {
        double xdiff, ydiff, norm;
        int shift;
        RealMatrix *X, *Y;
        
        if (j < n_top - 4) {
            shift = -2;
            X = &x_test_top;
            Y = &y_test_top;
        }
        else if (j < n_top + n_right - 8) {
            shift = n_top - 6;
            X = &x_test_right;
            Y = &y_test_right;
        }
        else if (j < n_top + n_right + n_bottom - 12) {
            shift = n_top + n_right - 10;
            X = &x_test_bottom;
            Y = &y_test_bottom;
        }
        else {
            shift = n_top + n_right + n_bottom - 14;
            X = &x_test_left;
            Y = &y_test_left;
        }
        xdiff = (*X).getHostValue(j - 1 - shift) - (*X).getHostValue(j + 1 - shift); // Central difference
        ydiff = (*Y).getHostValue(j - 1 - shift) - (*Y).getHostValue(j + 1 - shift); // Central difference
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;

        x_int.setHostValue(j, (*X).getHostValue(j - shift) - alpha*ydiff);
        y_int.setHostValue(j, (*Y).getHostValue(j - shift) + alpha*xdiff);

    }
    */
    
    bool save_segment = true;
    if (save_segment) {
        FILE *file;
        char filename[256];
        sprintf(filename,"../../../Data/segments/test_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_test; k++) {
            fprintf(file, "%.4e %.4e\n", x_test.getHostValue(k), y_test.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/ext_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_ext; k++) {
            fprintf(file, "%.4e %.4e\n", x_ext.getHostValue(k), y_ext.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/int_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_int; k++) {
            fprintf(file, "%.4e %.4e\n", x_int.getHostValue(k), y_int.getHostValue(k));
        }
        fclose(file);

    }


}





}