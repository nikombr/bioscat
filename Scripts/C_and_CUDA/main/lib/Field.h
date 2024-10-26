#ifndef _FIELD_H
#define _FIELD_H
extern "C" {
#include "ComplexMatrix.h"
struct Field {
    // Stores a real matrix rowmajor in a vector both on the host and the device
    public:
        ComplexMatrix x;
        ComplexMatrix y;
        ComplexMatrix z;

        // Methods
        Field() {
            x = ComplexMatrix();
            y = ComplexMatrix();
            z = ComplexMatrix();
        };  // Constructer
        Field(int rows) {
            x = ComplexMatrix(rows);
            y = ComplexMatrix(rows);
            z = ComplexMatrix(rows);
        };  // Constructer
        Field(int rows, bool x_bool, bool y_bool, bool z_bool) {
            if (x_bool) x = ComplexMatrix(rows);
            if (y_bool) y = ComplexMatrix(rows);
            if (z_bool) z = ComplexMatrix(rows);
        };  // Constructer
        Field(int rows, int cols) {
            x = ComplexMatrix(rows, cols);
            y = ComplexMatrix(rows, cols);
            z = ComplexMatrix(rows, cols);
        };  // Constructer
        Field(int rows, int cols, bool x_bool, bool y_bool, bool z_bool) {
            if (x_bool) x = ComplexMatrix(rows, cols);
            if (y_bool) y = ComplexMatrix(rows, cols);
            if (z_bool) z = ComplexMatrix(rows, cols);
        }; // Constructer
    

};
}

#endif