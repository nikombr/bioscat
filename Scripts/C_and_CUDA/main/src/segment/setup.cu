#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

void Segment::setup(Nanostructure nanostructure, int current_segment, int total_grid_points, int num_segments) {

    bool save_segment = false;

    int segment_length = total_grid_points / num_segments;
    int start, end, startnum, endnum;
    double startvalue, endvalue, step, startstep, endstep, startxvalue, endxvalue;
    int minNumSteps = 10;

    step = nanostructure.x.getHostValue(1) - nanostructure.x.getHostValue(0);
    double alpha = 2*step;

    FILE *file;
    char filename[256];

    start = current_segment * segment_length;
    end = min(start + segment_length + 1, total_grid_points);

    startvalue  = nanostructure.f.getHostValue(start);
    endvalue    = nanostructure.f.getHostValue(end - 1);

    startnum = max(minNumSteps, (int) ceil(startvalue/step));
    endnum   = max(minNumSteps, (int) ceil(endvalue/step));

    startstep = startvalue/startnum;
    endstep   = endvalue/endnum;
    
    startxvalue  = nanostructure.x.getHostValue(start);
    endxvalue    = nanostructure.x.getHostValue(end - 1);
    printf("HEJ!!\n");
    // Allocate arrays
    n_top      = end - start - 2;
    printf("ntop = %d",n_top);
    n_bottom   = end - start - 2;
    n_right    = endnum - 1;
    n_left     = startnum - 1;
    n_int = n_top + n_right + n_bottom + n_left - 16;
    n_ext = 2*(end - start - 1) + endnum + startnum;
    allocate(n_top, n_right, n_bottom, n_left, n_int, n_ext);
    printf("HEJ!!\n");
    RealMatrix x_temp = RealMatrix(n_ext + 2);
    RealMatrix y_temp = RealMatrix(n_ext + 2);

    // Save values
    num_test_points = n_top + n_bottom + n_right + n_left;
    num_interior_points = n_int;
    num_exterior_points = n_ext;

    

    // Compute temporary test points
    for (int j = 0; j < endnum; j++) {
        x_temp.setHostValue(end - 1 - start + j + 1, endxvalue);
        y_temp.setHostValue(end - 1 - start + n_right - j + 1, (j+1)*endstep);
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

    // Compute exterior points (eventuel lav bedre senere)
    //alpha = std::min(10e-8,2*std::min(startstep,endstep));
    for (int j = 1; j < n_ext + 1; j++) {
        double xdiff, ydiff, norm;

        xdiff = x_temp.getHostValue(j - 1) - x_temp.getHostValue(j + 1); // Central difference
        ydiff = y_temp.getHostValue(j - 1) - y_temp.getHostValue(j + 1); // Central difference
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;

        x_ext.setHostValue(j - 1, x_temp.getHostValue(j) + alpha*ydiff);
        y_ext.setHostValue(j - 1, y_temp.getHostValue(j) - alpha*xdiff);

        if (j >= 2 && j < n_top + 2) {
            n_x.setHostValue(j - 2,   ydiff);
            n_y.setHostValue(j - 2, - xdiff);
        }
    }

    // Remove end points
    start     += 1;
    end       -= 1;

    // Compute test points
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

 
    // Compute interior points
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

    
    if (save_segment) {
        sprintf(filename,"../../../Data/segments/test_top_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_top; k++) {
            fprintf(file, "%.4e %.4e\n", x_test_top.getHostValue(k), y_test_top.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_right_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_right; k++) {
            fprintf(file, "%.4e %.4e\n", x_test_right.getHostValue(k), y_test_right.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_bottom_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_bottom; k++) {
            fprintf(file, "%.4e %.4e\n", x_test_bottom.getHostValue(k), y_test_bottom.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_left_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_left; k++) {
            fprintf(file, "%.4e %.4e\n", x_test_left.getHostValue(k), y_test_left.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_int_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_int; k++) {
            fprintf(file, "%.4e %.4e\n", x_int.getHostValue(k), y_int.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_ext_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_ext; k++) {
            fprintf(file, "%.4e %.4e\n", x_ext.getHostValue(k), y_ext.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_temp_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_ext + 2; k++) {
            fprintf(file, "%.4e %.4e\n", x_temp.getHostValue(k), y_temp.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"../../../Data/segments/test_n_segment_%d.txt", current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_top; k++) {
            fprintf(file, "%.6f %.6f\n", n_x.getHostValue(k), n_y.getHostValue(k));
        }
        fclose(file);
    }

    x_temp.free();
    y_temp.free();

}


}