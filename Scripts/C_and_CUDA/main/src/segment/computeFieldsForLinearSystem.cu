#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

void Segment::computeFieldsForLinearSystem() {
    RealMatrix y = RealMatrix(num_test_points);
    RealMatrix x = RealMatrix(num_test_points);
    int shift = 0;
    for (int j = 0; j < y_test_top.rows;    j++) {
        x.setHostValue(j + shift,x_test_top.getHostValue(j));
        y.setHostValue(j + shift,y_test_top.getHostValue(j));
    }
    shift += y_test_top.rows;
    for (int j = 0; j < y_test_right.rows;  j++) {
        x.setHostValue(j + shift,x_test_right.getHostValue(j));
        y.setHostValue(j + shift,y_test_right.getHostValue(j));
    }
    shift += y_test_right.rows;
    for (int j = 0; j < y_test_bottom.rows; j++) {
        x.setHostValue(j + shift,x_test_bottom.getHostValue(j));
        y.setHostValue(j + shift,y_test_bottom.getHostValue(j));
    }
    shift += y_test_bottom.rows;
    for (int j = 0; j < y_test_left.rows;   j++) {
        x.setHostValue(j + shift,x_test_left.getHostValue(j));
        y.setHostValue(j + shift,y_test_left.getHostValue(j));
    }

    computeIncidentFieldVectors(y);
    computeReflectedFieldVectors(y);
    computeScatteredFieldMatrices(x, y, false);
    computeInteriorFieldMatrices(x, y);

    x.free();
    y.free();

}

}