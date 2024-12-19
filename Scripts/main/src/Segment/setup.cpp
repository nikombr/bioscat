#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#include "../../lib/Coordinates.h"
extern "C" {
using namespace std;

void computeExteriorPointsAndNormalVectors(Coordinates aux_ext, Coordinates normal_vectors, Coordinates aux_ext_temp, Nanostructure nanostructure, int start, int end, double alpha, double leftStep, double rightStep, int leftNum, int rightNum, int n_top, int n_right, int n_bottom, int n_left, double left_x_value, double right_x_value, int n_ext, bool deviceComputation, bool printOutput);
void computeTestPoints(Coordinates test_points, Nanostructure nanostructure, int start, int end, double leftStep, double rightStep,int n_top, int n_right, int n_bottom, int n_left, double left_x_value, double right_x_value, bool deviceComputation, bool printOutput);
void computeInteriorPoints(Coordinates aux_int,Coordinates test_points,double alpha, int n_top, int n_right, int n_bottom, int n_left, int n_int, bool deviceComputation, bool printOutput);

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
    computeExteriorPointsAndNormalVectors(aux_ext, normal_vectors, aux_ext_temp, nanostructure, start, end, alpha, leftStep, rightStep, leftNum, rightNum, n_top, n_right, n_bottom, n_left, left_x_value, right_x_value, n_ext, deviceComputation, printOutput);

    // Compute test points and temporary test points for interior points
    computeTestPoints(test_points, nanostructure, start,  end,  leftStep,  rightStep, n_top, n_right, n_bottom, n_left, left_x_value, right_x_value, deviceComputation, printOutput);

    // Compute interior points
    computeInteriorPoints(aux_int, test_points, alpha, n_top, n_right,  n_bottom,  n_left, n_int, deviceComputation, printOutput); 
    
    
    bool save_segment = false;
    if (save_segment) {
        /*test_points.allocateHost();
        aux_ext.allocateHost();
        aux_int.allocateHost();
        normal_vectors.allocateHost();*/
        if (deviceComputation) {
            test_points.toHost();
            aux_ext.toHost();
            aux_int.toHost();
            normal_vectors.toHost();
        }
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