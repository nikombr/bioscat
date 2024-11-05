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
        bool x_bool = false;
        bool y_bool = false;
        bool z_bool = false;

        // Methods
        Field() {
            x = ComplexMatrix();
            y = ComplexMatrix();
            z = ComplexMatrix();
        };  // Constructer
        Field(int rows) {
            x_bool = true;
            y_bool = true;
            z_bool = true;
            x = ComplexMatrix(rows);
            y = ComplexMatrix(rows);
            z = ComplexMatrix(rows);
        };  // Constructer
        Field(int rows, bool x_bool_new, bool y_bool_new, bool z_bool_new) {
            x_bool = x_bool_new;
            y_bool = y_bool_new;
            z_bool = z_bool_new;
            if (x_bool) x = ComplexMatrix(rows);
            if (y_bool) y = ComplexMatrix(rows);
            if (z_bool) z = ComplexMatrix(rows);
        };  // Constructer
        Field(int rows, int cols) {
            x_bool = true;
            y_bool = true;
            z_bool = true;
            x = ComplexMatrix(rows, cols);
            y = ComplexMatrix(rows, cols);
            z = ComplexMatrix(rows, cols);
        };  // Constructer
        Field(int rows, int cols, bool x_bool_new, bool y_bool_new, bool z_bool_new) {
            x_bool = x_bool_new;
            y_bool = y_bool_new;
            z_bool = z_bool_new;
            if (x_bool) x = ComplexMatrix(rows, cols);
            if (y_bool) y = ComplexMatrix(rows, cols);
            if (z_bool) z = ComplexMatrix(rows, cols);
        }; // Constructer

        void free() {
            if (x_bool) x.free();
            if (y_bool) y.free();
            if (z_bool) z.free();
        }

        void setHostZero() {
            if (x_bool) x.setHostZero();
            if (y_bool) y.setHostZero();
            if (z_bool) z.setHostZero();
        }

        void setDeviceZero() {
            if (x_bool) x.setDeviceZero();
            if (y_bool) y.setDeviceZero();
            if (z_bool) z.setDeviceZero();
        }
    

};
}

#endif