#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;


void Segment::setupRightHandSide() {
    // Setup right-hand side in linear system
    RealMatrix b_imag = RealMatrix(2 * num_test_points);
    RealMatrix b_real = RealMatrix(2 * num_test_points);

    double val, shift;
    //int imagShift = 2*num_test_points;
    ComplexMatrix *firstField_inc, * firstField_ref, *secondField_inc,*secondField_ref;
    //printf("scenario = %d\n",scenario);
    //printf("%d %\n",scenario == 1, scenario == 2);
    if (scenario == 1) {
        //printf("SCENARIO 1");
        firstField_inc = &E_inc_vector.z;
        firstField_ref = &E_ref_vector.z;
        secondField_inc = &H_inc_vector.x;
        secondField_ref = &H_ref_vector.x;
    }
    else if (scenario == 2) {
        firstField_inc = &H_inc_vector.x;
        firstField_ref = &H_ref_vector.x;
        secondField_inc = &E_inc_vector.z;
        secondField_ref = &E_ref_vector.z;
    }
    else {
        printf("You have to choose either 1 or 2 as the scenario!\n");
        return;
    }

    for (int j = 0; j < num_test_points; j++) {
        val =  - firstField_inc->getHostRealValue(j) - firstField_ref->getHostRealValue(j);
        b_real.setHostValue(j, val);
        val =  - firstField_inc->getHostImagValue(j) - firstField_ref->getHostImagValue(j);
        b_imag.setHostValue(j, val);
    }

    shift = num_test_points;
    for (int j = 0; j < n_top; j++) {
        val  = secondField_inc->getHostRealValue(j) + secondField_ref->getHostRealValue(j);
        val *= n_y.getHostValue(j);
        b_real.setHostValue(j + shift, val);
        val  = secondField_inc->getHostImagValue(j) + secondField_ref->getHostImagValue(j);
        val *= n_y.getHostValue(j);
        b_imag.setHostValue(j + shift, val);
    }

    shift += n_top;
    for(int j = 0; j < n_right; j++) {
        val = 0.0;
        b_real.setHostValue(j + shift, val);
        b_imag.setHostValue(j + shift, val);
    }

    shift += n_right;
    for(int j = 0; j < n_bottom; j++) {
        val  = secondField_inc->getHostRealValue(j + n_top + n_right) + secondField_ref->getHostRealValue(j + n_top + n_right);
        b_real.setHostValue(j + shift, val);
        val  = secondField_inc->getHostImagValue(j + n_top + n_right) + secondField_ref->getHostImagValue(j + n_top + n_right);
        val = 555;
        b_imag.setHostValue(j + shift, val);
    }

    shift += n_bottom;
    for(int j = 0; j < n_left; j++) {
        val = 0.0;
        b_real.setHostValue(j + shift, val);
        b_imag.setHostValue(j + shift, val);
    }
    /*printf("b:\n");
    for (int j = 0; j < 2*num_test_points; j++) {
        printf("%e\t + i(%e)\n",b_real.getHostValue(j),b_imag.getHostValue(j));
    }*/

    b = RealMatrix(4 * num_test_points);

    for (int r = 0; r < 2*num_test_points; r++) {
        b.setHostValue(r,                     b_real.getHostValue(r));
        b.setHostValue(r + 2*num_test_points, b_imag.getHostValue(r));
    }

    b_real.free();
    b_imag.free();
    
    return;
 
}

}