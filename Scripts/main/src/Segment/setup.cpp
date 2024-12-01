#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
extern "C" {
using namespace std;

void computeExteriorPointsAndNormalVectors(RealMatrix x_ext, RealMatrix y_ext, RealMatrix n_x, RealMatrix n_y, Nanostructure nanostructure, int start, int end, double alpha, double leftStep, double rightStep, int leftNum, int rightNum, int n_top, int n_right, int n_bottom, int n_left, double left_x_value, double right_x_value, bool deviceComputation, bool printOutput);
void computeTestPoints(RealMatrix x_test, RealMatrix y_test, Nanostructure nanostructure, int start, int end, double leftStep, double rightStep,int n_top, int n_right, int n_bottom, int n_left, double left_x_value, double right_x_value, bool deviceComputation, bool printOutput);
void computeInteriorPoints(RealMatrix x_int, RealMatrix y_int, RealMatrix x_test, RealMatrix y_test,double alpha, int n_top, int n_right, int n_bottom, int n_left, bool deviceComputation, bool printOutput);

void Segment::setup(Nanostructure nanostructure, int total_grid_points, int num_segments) {
    // Sets up test points and auxiliary points for the segment given a nanostructure
    // that is both on the host and on the device

    // Initialize variables
    int start, end, leftNum, rightNum;
    double leftStep, rightStep, left_f_value, right_f_value, left_x_value, right_x_value, alpha, step;

    // Determine distance of auxilliary points from nearest test point
    step = nanostructure.x.getHostValue(1) - nanostructure.x.getHostValue(0);

    // Determine where the segment starts and ends, we assume all segments have the same number of points
    start = current_segment * segment_length;
    end = start + segment_length;

    // Get values at end points
    left_f_value  = nanostructure.f.getHostValue(start);
    right_f_value = nanostructure.f.getHostValue(end - 1);
    left_x_value  = nanostructure.x.getHostValue(start);
    right_x_value = nanostructure.x.getHostValue(end - 1);

    // Determine steps for sides of segment
    leftNum   = n_left + 1;
    rightNum  = n_right + 1;
    leftStep  = left_f_value/leftNum;
    rightStep = right_f_value/rightNum;
    //alpha     = std::min((double)2*step,(double)1.0e-9),2*std::min(leftStep,rightStep);
    //alpha *=10;
    alpha = 1e-8;
    //printf("alpha = %e",alpha);
    
    // Compute exterior points
    computeExteriorPointsAndNormalVectors(aux_ext.x, aux_ext.y, normal_vectors.x, normal_vectors.y, nanostructure, start, end, alpha, leftStep, rightStep, leftNum, rightNum, n_top, n_right, n_bottom, n_left, left_x_value, right_x_value, deviceComputation, printOutput);

    // Compute test points and temporary test points for interior points
    computeTestPoints(test_points.x, test_points.y, nanostructure, start,  end,  leftStep,  rightStep, n_top, n_right, n_bottom, n_left, left_x_value, right_x_value, deviceComputation, printOutput);

    // Compute interior points
    computeInteriorPoints(aux_int.x, aux_int.y, test_points.x, test_points.y, alpha, n_top, n_right,  n_bottom,  n_left, deviceComputation, printOutput); 
    
    
    bool save_segment = true;
    if (save_segment) {
        test_points.allocateHost();
        aux_ext.allocateHost();
        aux_int.allocateHost();
        normal_vectors.allocateHost();
        test_points.toHost();
        aux_ext.toHost();
        aux_int.toHost();
        normal_vectors.toHost();
        FILE *file;
        char filename[256];
        char dir[256];
        sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");
        sprintf(filename,"%s/Data/segments/test_segment_%d.txt", dir, current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_test; k++) {
            fprintf(file, "%.4e %.4e\n", test_points.x.getHostValue(k), test_points.y.getHostValue(k));
        }
        fclose(file);
        
        sprintf(filename,"%s/Data/segments/ext_segment_%d.txt", dir, current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_ext; k++) {
            fprintf(file, "%.4e %.4e\n", aux_ext.x.getHostValue(k), aux_ext.y.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"%s/Data/segments/int_segment_%d.txt", dir, current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_int; k++) {
            fprintf(file, "%.4e %.4e\n", aux_int.x.getHostValue(k), aux_int.y.getHostValue(k));
        }
        fclose(file);

        sprintf(filename,"%s/Data/segments/n_segment_%d.txt", dir, current_segment+1);
        file = fopen(filename, "w");
        if (file == NULL) {
            perror("Error opening file");
            return;
        }
        for (int k = 0; k < n_test; k++) {
            fprintf(file, "%.4e %.4e\n", normal_vectors.x.getHostValue(k), normal_vectors.y.getHostValue(k));
        }
        fclose(file);

    }


}





}