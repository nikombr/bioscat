
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

    x_int.free();           y_int.free();
    x_ext.free();           y_ext.free();
    x_test_top.free();      y_test_top.free();
    x_test_right.free();    y_test_right.free();
    x_test_bottom.free();   y_test_bottom.free();
    x_test_left.free();     y_test_left.free();
    n_x.free();             n_y.free();

}

void Segment::allocate(int n_top, int n_right, int n_bottom, int n_left, int n_int, int n_ext) {

    x_int = RealMatrix(n_int);              y_int = RealMatrix(n_int);
    x_ext = RealMatrix(n_ext);              y_ext = RealMatrix(n_ext);
    x_test_top = RealMatrix(n_top);         y_test_top = RealMatrix(n_top);
    n_x = RealMatrix(n_top);                n_y = RealMatrix(n_top);
    x_test_right = RealMatrix(n_right);     y_test_right = RealMatrix(n_right);
    x_test_bottom = RealMatrix(n_bottom);   y_test_bottom = RealMatrix(n_bottom);
    x_test_left = RealMatrix(n_left);       y_test_left = RealMatrix(n_left);

}


void Segment::setup(Nanostructure nanostructure, int current_segment, int total_grid_points, int num_segments) {

    bool save_segment = true;

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

    // Save values
    num_test_points = n_top + n_bottom + n_right + n_left;
    num_interior_points = n_int;
    num_exterior_points = n_ext;

    // Compute test points
    for (int j = 0; j < endnum; j++) {
        x_temp.setHostValue(end - 1 - start + j, endxvalue);
        y_temp.setHostValue(end - 1 - start + n_right - j, (j+1)*endstep);
    }
    
    for (int j = 0; j < startnum + 1; j++) {
        x_temp.setHostValue(2*(end - start - 1) + endnum + j, startxvalue);
        y_temp.setHostValue(2*(end - start - 1) + endnum + j, (j+1)*startstep);
    }

    for (int j = start; j < end - 1; j++) {
        
        x_temp.setHostValue(j - start, nanostructure.x.getHostValue(j));
        x_temp.setHostValue(end - start - 1 + endnum + end - j - 2, nanostructure.x.getHostValue(j));
        y_temp.setHostValue(j - start, nanostructure.f.getHostValue(j));
        y_temp.setHostValue(end - start - 1 + endnum + j - start + 1 - 1, 0.0);
    }

    // Compute exterior points (eventuel lav bedre senere)
    for (int j = 0; j < n_ext; j++) {
        double xdiff, ydiff, norm;

        xdiff = x_temp.getHostValue(j) - x_temp.getHostValue(j + 1);
        ydiff = y_temp.getHostValue(j) - y_temp.getHostValue(j + 1);
        
        norm = std::sqrt(xdiff*xdiff + ydiff*ydiff);
        xdiff /= norm;
        ydiff /= norm;

        x_ext.setHostValue(j, x_temp.getHostValue(j) + alpha*ydiff);
        y_ext.setHostValue(j, y_temp.getHostValue(j) - alpha*xdiff);

        if (j > 0 && j < n_top + 1) {
            n_x.setHostValue(j-1,ydiff);
            n_y.setHostValue(j-1,-xdiff);
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
        xdiff = (*X).getHostValue(j - shift) - (*X).getHostValue(j + 1 - shift);
        ydiff = (*Y).getHostValue(j - shift) - (*Y).getHostValue(j + 1 - shift);
        
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
        for (int k = 0; k < n_ext; k++) {
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
            fprintf(file, "%.4e %.4e\n", n_x.getHostValue(k), n_y.getHostValue(k));
        }
        fclose(file);
    }

    x_temp.free();
    y_temp.free();

}

double H02_real(double x) {
    int n = 0;
    // Compute Bessel functions of the first (Jn) and second (Yn) kinds
    double Jn = jn(n, x);
    // Hankel functions
    return Jn;          // Real part of H_n^(2)(x)
}

double H02_imaginary(double x) {
    int n = 0;
    // Compute Bessel functions of the first (Jn) and second (Yn) kinds
    double Yn = yn(n, x);
    // Hankel functions
    return -Yn;         // Imaginary part of H_n^(2)(x)
}

double H12_real(double x) {
    int n = 1;
    // Compute Bessel functions of the first (Jn) and second (Yn) kinds
    double Jn = jn(n, x);
    // Hankel functions
    return Jn;          // Real part of H_n^(2)(x)
}

double H12_imaginary(double x) {
    int n = 1;
    // Compute Bessel functions of the first (Jn) and second (Yn) kinds
    double Yn = yn(n, x);
    // Hankel functions
    return -Yn;         // Imaginary part of H_n^(2)(x)
}

void Segment::computeScatteredFieldMatrices(RealMatrix x, RealMatrix y, bool far_field_approximation) {
    
    int rows = y.rows;
    int cols = y_int.rows;
    double abs_int, abs_int_ref, xdiff, ydiff, ydiff_ref, H_real, H_imaginary, H_real_ref, H_imaginary_ref, val;

    // Determine which matrices to be allocated
    bool Ex_bool = scenario == 1 ? false : true;
    bool Ey_bool = scenario == 1 ? false : true;
    bool Ez_bool = scenario == 1 ? true  : false;
    bool Hx_bool = scenario == 1 ? true  : false;
    bool Hy_bool = scenario == 1 ? true  : false;
    bool Hz_bool = scenario == 1 ? false : true;

    E_scat_matrix = Field(rows, cols, Ex_bool, Ey_bool, Ez_bool);
    H_scat_matrix = Field(rows, cols, Hx_bool, Hy_bool, Hz_bool);

    if (scenario == 1) {
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {

                // Get data
                xdiff       = x.getHostValue(r) - x_int.getHostValue(c);
                ydiff       = y.getHostValue(r) - y_int.getHostValue(c);
                ydiff_ref   = y.getHostValue(r) + y_int.getHostValue(c);
                abs_int     = std::sqrt(xdiff*xdiff + ydiff*ydiff);
                abs_int_ref = std::sqrt(xdiff*xdiff + ydiff_ref*ydiff_ref);

                // Compute first Hankel functions
                H_real          = H02_real(constants.k0*abs_int);
                H_real_ref      = H02_real(constants.k0*abs_int_ref);
                H_imaginary     = H02_imaginary(constants.k0*abs_int);
                H_imaginary_ref = H02_imaginary(constants.k0*abs_int_ref);
                
                val = H_real + constants.Gamma_ref * H_real_ref;
                E_scat_matrix.z.setHostRealValue(r, c, val);
                val = H_imaginary + constants.Gamma_ref * H_imaginary_ref;
                E_scat_matrix.z.setHostImaginaryValue(r, c, val);

                // Compute second Hankel functions
                H_real          = H12_real(constants.k0*abs_int);
                H_real_ref      = H12_real(constants.k0*abs_int_ref);
                H_imaginary     = H12_imaginary(constants.k0*abs_int);
                H_imaginary_ref = H12_imaginary(constants.k0*abs_int_ref);

                val = 1/constants.eta0 * (1/abs_int      * H_imaginary     * ydiff + \
                     constants.Gamma_ref * 1/abs_int_ref * H_imaginary_ref * ydiff_ref);
                H_scat_matrix.x.setHostRealValue(r, c, val);
                val = -1/constants.eta0 * (1/abs_int     * H_real     * ydiff + \
                     constants.Gamma_ref * 1/abs_int_ref * H_real_ref * ydiff_ref);
                H_scat_matrix.x.setHostImaginaryValue(r, c, val);

                val = -1/constants.eta0 * xdiff * (1/abs_int      * H_imaginary      + \
                             constants.Gamma_ref * 1/abs_int_ref  * H_imaginary_ref);
                H_scat_matrix.x.setHostRealValue(r, c, val);
                val = 1/constants.eta0 * xdiff * (1/abs_int     * H_real      + \
                            constants.Gamma_ref * 1/abs_int_ref * H_real_ref);
                H_scat_matrix.x.setHostImaginaryValue(r, c, val);
            }
        }

    }
    else if (scenario == 2) {

    }
    else {
        printf("Please input 1 or 2 for the scenario!\n");
    }

}

void Segment::computeInteriorFieldMatrices(RealMatrix x, RealMatrix y) {
    
    int rows = y.rows;
    int cols = y_ext.rows;
    double abs_ext, xdiff, ydiff, H_real, H_imaginary, val;

    // Determine which matrices to be allocated
    bool Ex_bool = scenario == 1 ? false : true;
    bool Ey_bool = scenario == 1 ? false : true;
    bool Ez_bool = scenario == 1 ? true  : false;
    bool Hx_bool = scenario == 1 ? true  : false;
    bool Hy_bool = scenario == 1 ? true  : false;
    bool Hz_bool = scenario == 1 ? false : true;

    E_int_matrix = Field(rows, cols, Ex_bool, Ey_bool, Ez_bool);
    H_int_matrix = Field(rows, cols, Hx_bool, Hy_bool, Hz_bool);

    if (scenario == 1) {
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {

                // Get data
                xdiff       = x.getHostValue(r) - x_ext.getHostValue(c);
                ydiff       = y.getHostValue(r) - y_ext.getHostValue(c);
                abs_ext     = std::sqrt(xdiff*xdiff + ydiff*ydiff);

                // Compute first Hankel functions
                H_real      = H02_real(constants.k1*abs_ext);
                H_imaginary = H02_imaginary(constants.k1*abs_ext);
                
                val = H_real;
                E_scat_matrix.z.setHostRealValue(r, c, val);
                val = H_imaginary;
                E_scat_matrix.z.setHostImaginaryValue(r, c, val);

                // Compute second Hankel functions
                H_real          = H12_real(constants.k0*abs_ext);
                H_imaginary     = H12_imaginary(constants.k0*abs_ext);

                val = constants.n1/constants.eta0*1/abs_ext*ydiff*H_imaginary;
                H_scat_matrix.x.setHostRealValue(r, c, val);
                val = -constants.n1/constants.eta0*1/abs_ext*ydiff*H_real;
                H_scat_matrix.x.setHostImaginaryValue(r, c, val);

                val = -constants.n1/constants.eta0*1/abs_ext*xdiff*H_imaginary;
                H_scat_matrix.x.setHostRealValue(r, c, val);
                val = constants.n1/constants.eta0*1/abs_ext*xdiff*H_imaginary;
                H_scat_matrix.x.setHostImaginaryValue(r, c, val);
            }
        }

    }
    else if (scenario == 2) {

    }
    else {
        printf("Please input 1 or 2 for the scenario!\n");
    }

}

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
        x.setHostValue(j + shift,x_test_top.getHostValue(j));
        y.setHostValue(j + shift,y_test_top.getHostValue(j));
    }
    shift += y_test_right.rows;
    for (int j = 0; j < y_test_bottom.rows; j++) {
        x.setHostValue(j + shift,x_test_top.getHostValue(j));
        y.setHostValue(j + shift,y_test_top.getHostValue(j));
    }
    shift += y_test_bottom.rows;
    for (int j = 0; j < y_test_left.rows;   j++) {
        x.setHostValue(j + shift,x_test_top.getHostValue(j));
        y.setHostValue(j + shift,y_test_top.getHostValue(j));
    }

    computeIncidentFieldVectors(y);
    computeReflectedFieldVectors(y);
    computeScatteredFieldMatrices(x, y, false);
    computeInteriorFieldMatrices(x, y);

    y.free();

}


void Segment::setupSystemSubMatrices() {

    A_real = RealMatrix(2 * num_test_points, n_int + n_ext);
    A_imaginary = RealMatrix(2 * num_test_points, n_int + n_ext);

    double val, rshift;
    ComplexMatrix * firstField_scat,  * firstField_int,  \
                  * secondField_scat, * secondField_int, \
                  * thirdField_scat,  * thirdField_int;

    if (scenario == 1) {
        firstField_scat  = &E_scat_matrix.z;
        firstField_int   = &E_scat_matrix.z;
        secondField_scat = &H_scat_matrix.x;
        secondField_int  = &H_scat_matrix.x;
        thirdField_scat  = &H_scat_matrix.y;
        thirdField_int   = &H_scat_matrix.y;
    }
    else if (scenario == 2) {
        firstField_scat  = &H_scat_matrix.z;
        firstField_int   = &H_scat_matrix.z;
        secondField_scat = &E_scat_matrix.x;
        secondField_int  = &E_scat_matrix.x;
        thirdField_scat  = &E_scat_matrix.y;
        thirdField_int   = &E_scat_matrix.y;
    }
    else {
        printf("You have to choose either 1 or 2 as the scenario!\n");
        return;
    }

    for (int r = 0; r < num_test_points; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  firstField_scat->getHostRealValue(r, c);
            A_real.setHostValue(r, c, val);
            val =  firstField_scat->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -firstField_int->getHostRealValue(r, c);
            A_real.setHostValue(r, c + n_int, val);
            val =  -firstField_int->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r, c + n_int, val);
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
            val += - n_y.getHostValue(r) * secondField_scat->getHostImaginaryValue(r, c);
            val +=   n_x.getHostValue(r) * thirdField_scat->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val  = 0;
            val +=   n_y.getHostValue(r) * secondField_int->getHostRealValue(r, c);
            val += - n_x.getHostValue(r) * thirdField_int->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val  = 0;
            val +=   n_y.getHostValue(r) * secondField_int->getHostImaginaryValue(r, c);
            val += - n_x.getHostValue(r) * thirdField_int->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c + n_int, val);
        }
    }

    rshift += n_top;
    for(int r = 0; r < n_right; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  thirdField_scat->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c, val);
            val =  thirdField_scat->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -thirdField_int->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val =  -thirdField_int->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c + n_int, val);
        }
    }

    rshift += n_right;
    for(int r = 0; r < n_bottom; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  secondField_scat->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c, val);
            val =  secondField_scat->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -secondField_int->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val =  -secondField_int->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c + n_int, val);
        }
    }

    rshift += n_bottom;
    for(int r = 0; r < n_left; r++) {
        for (int c = 0; c < n_int; c++) {
            val =  thirdField_scat->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c, val);
            val =  thirdField_scat->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c, val);
        }
        for (int c = 0; c < n_ext; c++) {
            val =  -thirdField_int->getHostRealValue(r, c);
            A_real.setHostValue(r + rshift, c + n_int, val);
            val =  -thirdField_int->getHostImaginaryValue(r, c);
            A_imaginary.setHostValue(r + rshift, c + n_int, val);
        }
    }
    
    return;
 
}


}