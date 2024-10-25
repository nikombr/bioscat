
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/2D/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;


Segment::Segment() {
    // Empty constructor
}
    
void Segment::free() {
    x_int.free();
    y_int.free();
    x_ext.free();
    y_ext.free();
    x_test_top.free();
    y_test_top.free();
    x_test_right.free();
    y_test_right.free();
    x_test_bottom.free();
    y_test_bottom.free();
    x_test_left.free();
    y_test_left.free();
}

void Segment::allocate(int n_top, int n_right, int n_bottom, int n_left, int n_int, int n_ext) {

    x_int = RealMatrix(n_int);
    y_int = RealMatrix(n_int);
    x_ext = RealMatrix(n_ext);
    y_ext = RealMatrix(n_ext);
    x_test_top = RealMatrix(n_top);
    y_test_top = RealMatrix(n_top);
    x_test_right = RealMatrix(n_right);
    y_test_right = RealMatrix(n_right);
    x_test_bottom = RealMatrix(n_bottom);
    y_test_bottom = RealMatrix(n_bottom);
    x_test_left = RealMatrix(n_left);
    y_test_left = RealMatrix(n_left);

}


void Segment::setup(Nanostructure nanostructure, int current_segment, int total_grid_points, int num_segments) {

    bool save_segment = true;

    int segment_length = total_grid_points / num_segments;
    int start, end, startnum, endnum;
    double startvalue, endvalue, step, startstep, endstep, startxvalue, endxvalue;
    int n_top,  n_right,  n_bottom,  n_left,  n_int,  n_ext;
    int minNumSteps = 10;

    step = nanostructure.x.getHostValue(1) - nanostructure.x.getHostValue(0);
    double alpha = 2*step;

    FILE *file;
    char filename[256];
    //printf("hej %d\n",current_segment);

    start = current_segment * segment_length;
    end = min(start + segment_length + 1, total_grid_points);


    startvalue  = nanostructure.f.getHostValue(start);
    endvalue    = nanostructure.f.getHostValue(end - 1);

    startnum = max(minNumSteps, (int) ceil(startvalue/step));
    endnum   = max(minNumSteps, (int) ceil(endvalue/step));
    //printf("hmm = %f\n",ceil(startvalue/step));

    startstep = startvalue/startnum;
    endstep   = endvalue/endnum;

    startxvalue  = nanostructure.x.getHostValue(start);
    endxvalue    = nanostructure.x.getHostValue(end - 1);

    // Allocate arrays
    n_top      = end - start - 2;
    n_bottom   = end - start - 2;
    n_right    = endnum - 1;
    n_left     = startnum - 1;
    n_int = n_top + n_right + n_bottom + n_left - 16;
    n_ext = 2*(end - start - 1) + endnum + startnum;
    allocate(n_top, n_right, n_bottom, n_left, n_int, n_ext);
    RealMatrix x_temp = RealMatrix(n_ext + 1);
    RealMatrix y_temp = RealMatrix(n_ext + 1);

    // Compute test points
    for (int j = 0; j < endnum; j++) {
        x_temp.setHostValue(end - start + j, endxvalue);
        y_temp.setHostValue(end - start + n_right - j - 1, (j+1)*endstep);
    }
    
    for (int j = 0; j < startnum; j++) {
        x_temp.setHostValue(2*(end - start - 1) + endnum + j, startxvalue);
        y_temp.setHostValue(2*(end - start - 1) + endnum + j, j*startstep);
    }

    x_temp.setHostValue(n_ext, startxvalue);
    y_temp.setHostValue(n_ext, startnum*startstep);

    for (int j = start; j < end - 1; j++) {
        
        x_temp.setHostValue(j - start, nanostructure.x.getHostValue(j));
        y_temp.setHostValue(j - start, nanostructure.f.getHostValue(j));
        x_temp.setHostValue(end - start - 1 + endnum + end - j, nanostructure.x.getHostValue(j));
        y_temp.setHostValue(end - start - 1 + endnum + j - start + 1, 0.0);
    }

    // Compute exterior points
    for (int j = 0; j < n_ext; j++) {
        double xdiff, ydiff, norm;
        int shift;
        RealMatrix *X, *Y;
        
        if (j < n_top - 1) {
            shift = 0;
            X = &x_test_top;
            Y = &y_test_top;
        }
        else if (j < n_top + n_right - 2) {
            shift = n_top - 1;
            X = &x_test_right;
            Y = &y_test_right;
        }
        else if (j < n_top + n_right + n_bottom - 3) {
            shift = n_top + n_right - 2;
            X = &x_test_bottom;
            Y = &y_test_bottom;
        }
        else {
            shift = n_top + n_right + n_bottom - 3;
            X = &x_test_left;
            Y = &y_test_left;
        }
        printf("hej %d %f\n",j-shift, (*X).getHostValue(j + 1 - shift));

        xdiff = (*X).getHostValue(j - shift) - (*X).getHostValue(j + 1 - shift);
        ydiff = (*Y).getHostValue(j - shift) - (*Y).getHostValue(j + 1 - shift);

        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff *= alpha/norm;
        ydiff *= alpha/norm;

        x_ext.setHostValue(j, (*X).getHostValue(j - shift) + ydiff);
        y_ext.setHostValue(j, (*Y).getHostValue(j - shift) - xdiff);
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
        
        if (j < n_top - 4) {
            shift = -2;
            xdiff = x_test_top.getHostValue(j - shift) - x_test_top.getHostValue(j + 1 - shift);
            ydiff = y_test_top.getHostValue(j - shift) - y_test_top.getHostValue(j + 1 - shift);
        }
        else if (j < n_top + n_right - 8) {
            shift = n_top - 6;
            xdiff = x_test_right.getHostValue(j - shift) - x_test_right.getHostValue(j + 1 - shift);
            ydiff = y_test_right.getHostValue(j - shift) - y_test_right.getHostValue(j + 1 - shift);
        }
        else if (j < n_top + n_right + n_bottom - 12) {
            shift = n_top + n_right - 10;
            xdiff = x_test_bottom.getHostValue(j - shift) - x_test_bottom.getHostValue(j + 1 - shift);
            ydiff = y_test_bottom.getHostValue(j - shift) - y_test_bottom.getHostValue(j + 1 - shift);
        }
        else {
            shift = n_top + n_right + n_bottom - 14;
            xdiff = x_test_left.getHostValue(j - shift) - x_test_left.getHostValue(j + 1 - shift);
            ydiff = y_test_left.getHostValue(j - shift) - y_test_left.getHostValue(j + 1 - shift);
        }
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff *= alpha/norm;
        ydiff *= alpha/norm;

        if (j < n_top - 4) {
            shift = -2;
            x_int.setHostValue(j, x_test_top.getHostValue(j - shift) - ydiff);
            y_int.setHostValue(j, y_test_top.getHostValue(j - shift) + xdiff);
        }
        else if (j < n_top + n_right - 8) {
            shift = n_top - 6;
            x_int.setHostValue(j, x_test_right.getHostValue(j - shift) - ydiff);
            y_int.setHostValue(j, y_test_right.getHostValue(j - shift) + xdiff);
        }
        else if (j < n_top + n_right + n_bottom - 12) {
            shift = n_top + n_right - 10;
            x_int.setHostValue(j, x_test_bottom.getHostValue(j - shift) - ydiff);
            y_int.setHostValue(j, y_test_bottom.getHostValue(j - shift) + xdiff);
        }
        else {
            shift = n_top + n_right + n_bottom - 14;
            x_int.setHostValue(j, x_test_left.getHostValue(j - shift) - ydiff);
            y_int.setHostValue(j, y_test_left.getHostValue(j - shift) + xdiff);
        }
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
        for (int k = 0; k < n_ext; k++) {
            fprintf(file, "%.4e %.4e\n", x_temp.getHostValue(k), y_temp.getHostValue(k));
        }
        fclose(file);
    }

    x_temp.free();
    y_temp.free();

}




}