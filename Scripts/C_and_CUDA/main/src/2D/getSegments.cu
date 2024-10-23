#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/2D/Segment.h"
#include "../../lib/BioScat.h"
using namespace std;


void BioScat::getSegments(int num_segments) {

    this->segments = new Segment[num_segments];

    int segment_length = total_grid_points / num_segments;
    int start, end, startnum, endnum;
    double startvalue, endvalue, step, startstep, endstep, startxvalue, endxvalue;
    int n1,  n2,  n3,  n4,  n_int,  n_ext;
    int minNumSteps = 10;

    step = nanostructure.f.val_h[1] - nanostructure.f.val_h[0];

    for (int i = 0; i < num_segments; i++) {

        start = i * segment_length;
        end = min(start + segment_length + 1, total_grid_points);

        startvalue  = nanostructure.f.val_h[start];
        endvalue    = nanostructure.f.val_h[end - 1];

        startnum = max(minNumSteps, (int) ceil(step/startvalue));
        endnum   = max(minNumSteps, (int) ceil(step/endvalue));

        startstep = startvalue/startnum;
        endstep   = endvalue/endnum;

        // Allocate arrays
        n1    = segment_length - 2;
        n3    = segment_length - 2;
        n2    = endnum - 2;
        n4    = startnum - 2;
        n_int = 2*segment_length + endnum + startnum - 8;
        n_ext = 2*segment_length + endnum + startnum - 8;
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
            
            segments[i].x_test_top.val_h[j - start] = nanostructure.f.val_h[j];
            segments[i].y_test_top.val_h[j - start] = nanostructure.x.val_h[j];
            segments[i].x_test_bottom.val_h[end - j - 1 - start] = nanostructure.f.val_h[j];
            segments[i].y_test_bottom.val_h[j - start] = 0;
        }

    }



}

}