#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;


void Segment::setupSystemMatrix() {

    double val, rshift;
    ComplexMatrix * firstField_scat,  * firstField_int,  \
                  * secondField_scat, * secondField_int, \
                  * thirdField_scat,  * thirdField_int;

    if (scenario == 1) {
        firstField_scat  = &E_scat_matrix.z;
        firstField_int   = &E_int_matrix.z;
        secondField_scat = &H_scat_matrix.x;
        secondField_int  = &H_int_matrix.x;
        thirdField_scat  = &H_scat_matrix.y;
        thirdField_int   = &H_int_matrix.y;
    }
    else if (scenario == 2) {
        firstField_scat  = &H_scat_matrix.z;
        firstField_int   = &H_int_matrix.z;
        secondField_scat = &E_scat_matrix.x;
        secondField_int  = &E_int_matrix.x;
        thirdField_scat  = &E_scat_matrix.y;
        thirdField_int   = &E_int_matrix.y;
    }
    else {
        printf("You have to choose either 1 or 2 as the scenario!\n");
        return;
    }

    int num = firstField_scat->rows; // Number of test or observation points

    RealMatrix A_real = RealMatrix(2 * num_test_points, n_int + n_ext);
    RealMatrix A_imag = RealMatrix(2 * num_test_points, n_int + n_ext);

    for (int r = 0; r < num_test_points; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  firstField_scat->getHostRealValue(r, c);
            A_real.setHostValue(r, c, val);
            val =  firstField_scat->getHostImagValue(r, c);
            A_imag.setHostValue(r, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -firstField_int->getHostRealValue(r, c);
            A_real.setHostValue(r, c + n_int, val);
            val =  -firstField_int->getHostImagValue(r, c);
            A_imag.setHostValue(r, c + n_int, val);
        }
       
    }

    rshift = num_test_points;
    for (int r = 0; r < n_top; r++) {
        for (int c = 0; c < n_int; c++) {
            val  = 0;
            val += - n_y.getHostValue(r) * secondField_scat->getHostRealValue(r, c);
            val +=   n_x.getHostValue(r) * thirdField_scat->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c, val);
            val  = 0;
            val += - n_y.getHostValue(r) * secondField_scat->getHostImagValue(r, c);
            val +=   n_x.getHostValue(r) * thirdField_scat->getHostImagValue(r, c);
            A_imag.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val  = 0;
            val +=   n_y.getHostValue(r) * secondField_int->getHostRealValue(r, c);
            val += - n_x.getHostValue(r) * thirdField_int->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val  = 0;
            val +=   n_y.getHostValue(r) * secondField_int->getHostImagValue(r, c);
            val += - n_x.getHostValue(r) * thirdField_int->getHostImagValue(r, c);
            A_imag.setHostValue(r + rshift, c + n_int, val);
        }
    }

    rshift += n_top;
    for(int r = 0; r < n_right; r++) {
        for (int c = 0; c < n_int; c++) {
            
            val =  thirdField_scat->getHostRealValue(r + n_top, c);
            //val = 33.3;
            A_real.setHostValue(r + rshift, c, val);
            val =  thirdField_scat->getHostImagValue(r + n_top, c);
            A_imag.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -thirdField_int->getHostRealValue(r + n_top, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val =  -thirdField_int->getHostImagValue(r + n_top, c);
            A_imag.setHostValue(r + rshift, c + n_int, val);
        }
    }

    rshift += n_right;
    for(int r = 0; r < n_bottom; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  secondField_scat->getHostRealValue(r + n_top + n_right, c);
            A_real.setHostValue(r + rshift, c, val);
            val =  secondField_scat->getHostImagValue(r + n_top + n_right, c);
            A_imag.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -secondField_int->getHostRealValue(r + n_top + n_right, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val =  -secondField_int->getHostImagValue(r + n_top + n_right, c);
            A_imag.setHostValue(r + rshift, c + n_int, val);
        }
    }

    rshift += n_bottom;
    for(int r = 0; r < n_left; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  thirdField_scat->getHostRealValue(r + n_top + n_right + n_bottom, c);
            A_real.setHostValue(r + rshift, c, val);
            val =  thirdField_scat->getHostImagValue(r + n_top + n_right + n_bottom, c);
            A_imag.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -thirdField_int->getHostRealValue(r + n_top + n_right + n_bottom, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val =  -thirdField_int->getHostImagValue(r + n_top + n_right + n_bottom, c);
            A_imag.setHostValue(r + rshift, c + n_int, val);
        }
    }

    /*FILE *file;
    char * filename = "A_real_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < 2*num_test_points; r++) {
        for (int c = 0; c < n_int + n_ext; c++) {
            fprintf(file, "%e ", A_real.getHostValue(r,c));
        }
        fprintf(file, "\n");
    }
    fclose(file);
    filename = "A_imag_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < 2*num_test_points; r++) {
        for (int c = 0; c < n_int + n_ext; c++) {
            fprintf(file, "%e ", A_imag.getHostValue(r,c));
        }
        fprintf(file, "\n");
    }
    fclose(file);


    filename = "b_real_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < 2*num_test_points; r++) {
        fprintf(file, "%e\n", b_real.getHostValue(r));
    }
    fclose(file);
    filename = "b_imag_C.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < 2*num_test_points; r++) {
        fprintf(file, "%e\n", b_imag.getHostValue(r));
    }
    fclose(file);*/
    
    A = RealMatrix(4 * num_test_points, 2*(n_ext + n_int));
    for (int r = 0; r < 2 * num_test_points; r++) {
        for (int c = 0; c < n_ext + n_int; c++) {
            A.setHostValue(r,                       c,                   A_real.getHostValue(r,c));
            A.setHostValue(r,                       c + n_ext + n_int, - A_imag.getHostValue(r,c));
            A.setHostValue(r + 2 * num_test_points, c,                   A_imag.getHostValue(r,c));
            A.setHostValue(r + 2 * num_test_points, c + n_ext + n_int,   A_real.getHostValue(r,c));
        }
    }

    A_real.free();
    A_imag.free();

    return;
}


}