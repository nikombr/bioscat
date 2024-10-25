#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/2D/Segment.h"
#include "../../lib/BioScat.h"
#include "../../lib/RealMatrix.h"
using namespace std;


void BioScat::getSegments() {

    bool save_segment = true;

    this->segments = new Segment[num_segments];

    int segment_length = total_grid_points / num_segments;
    int start, end, startnum, endnum;
    double startvalue, endvalue, step, startstep, endstep, startxvalue, endxvalue;
    int n1,  n2,  n3,  n4,  n_int,  n_ext;
    int minNumSteps = 10;

    step = nanostructure.x.val_h[1] - nanostructure.x.val_h[0];

    FILE *file;
    char filename[256];

    for (int i = 0; i < num_segments; i++) {
        printf("hej %d\n",i);

        start = i * segment_length;
        end = min(start + segment_length + 1, total_grid_points);


        startvalue  = nanostructure.f.val_h[start];
        endvalue    = nanostructure.f.val_h[end - 1];

        startnum = max(minNumSteps, (int) ceil(step/startvalue));
        endnum   = max(minNumSteps, (int) ceil(step/endvalue));

        startstep = startvalue/startnum;
        endstep   = endvalue/endnum;

        // Allocate arrays
        n1    = end - start - 2;
        n3    = end - start - 2;
        n2    = endnum - 2;
        n4    = startnum - 2;
        n_int = 2*(end - start) + endnum + startnum - 8;
        n_ext = 2*(end - start) + endnum + startnum - 8;
        segments[i] = Segment(n1, n2, n3, n4, n_int, n_ext);

        startxvalue  = nanostructure.x.val_h[start];
        endxvalue    = nanostructure.x.val_h[end - 1];

        // Remove end points
        start     += 1;
        end       -= 1;

        for (int j = 0; j < startnum - 2; j++) {
            segments[i].x_test_left.val_h[j]  = startxvalue;
            segments[i].x_test_right.val_h[j] = endxvalue;
            segments[i].y_test_left.val_h[j]  = (j+1)*startstep;
            segments[i].y_test_right.val_h[j] = (j+1)*endstep;
        }

        for (int j = start; j < end; j++) {
            
            segments[i].x_test_top.val_h[j - start] = nanostructure.x.val_h[j];
            segments[i].y_test_top.val_h[j - start] = nanostructure.f.val_h[j];
            segments[i].x_test_bottom.val_h[end - j - 1] = nanostructure.x.val_h[j];
            printf("%d %d %d\n",j - start,end - j - 1,j);
            segments[i].y_test_bottom.val_h[j - start] = 0;
        }
        
        if (save_segment) {
            sprintf(filename,"../../../Data/segments/test_top_segment_%d.txt", i+1);
            file = fopen(filename, "w");
            if (file == NULL) {
                perror("Error opening file");
                return;
            }
            for (int k = 0; k < n1; k++) {
                fprintf(file, "%.4e %.4e\n", segments[i].x_test_top.val_h[k], segments[i].y_test_top.val_h[k]);
            }
            fclose(file);

            sprintf(filename,"../../../Data/segments/test_right_segment_%d.txt", i+1);
            file = fopen(filename, "w");
            if (file == NULL) {
                perror("Error opening file");
                return;
            }
            for (int k = 0; k < n2; k++) {
                fprintf(file, "%.4e %.4e\n", segments[i].x_test_right.val_h[k], segments[i].y_test_right.val_h[k]);
            }
            fclose(file);

            sprintf(filename,"../../../Data/segments/test_bottom_segment_%d.txt", i+1);
            file = fopen(filename, "w");
            if (file == NULL) {
                perror("Error opening file");
                return;
            }
            for (int k = 0; k < n3; k++) {
                fprintf(file, "%.4e %.4e\n", segments[i].x_test_bottom.val_h[k], segments[i].y_test_bottom.val_h[k]);
            }
            fclose(file);

            sprintf(filename,"../../../Data/segments/test_left_segment_%d.txt", i+1);
            file = fopen(filename, "w");
            if (file == NULL) {
                perror("Error opening file");
                return;
            }
            for (int k = 0; k < n4; k++) {
                fprintf(file, "%.4e %.4e\n", segments[i].x_test_left.val_h[k], segments[i].y_test_left.val_h[k]);
            }
            fclose(file);
        }
    }



}

}