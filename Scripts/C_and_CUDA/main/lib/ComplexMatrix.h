#ifndef _COMPLEX_MATRIX_H
#define _COMPLEX_MATRIX_H

class ComplexMatrix {
    // Stores a imaginaryMatrix rowmajor in a vector both on the host and the device
    private: // Private to avoid wrong indexing
        double * real_h;    // Real entries on host
        double * imaginary_h; // Imaginary entries on host
        double * real_d;    // Real entries on device
        double * imaginary_d; // Imaginary entries on device
    public:
        int    rows; // Number of rows
        int    cols; // Number of cols

        // Methods
        ComplexMatrix() {
            real_h = NULL;
            imaginary_h = NULL;
            real_d = NULL;
            imaginary_d = NULL;
        } 
        ComplexMatrix(int rows, int cols);                                  // Constructer, allocates arrays and initializes to zero
        ComplexMatrix(int rows);
        void free();                                           // Frees arrays
        void toHost();                                                      // Sends data to host
        void toDevice();                                                    // Sends data to device
        void setHostRealValue(int r, double val);                           // Sets host real value for vectors
        void setHostRealValue(int r, int c, double val);                    // Sets host real value for matrices
        void setHostImaginaryValue(int r, double val);                        // Sets host imaginary value for vectors
        void setHostImaginaryValue(int r, int c, double val);                 // Sets host imaginary value for matrices
        __device__ void setDeviceRealValue(int r, double val);              // Sets device real value for vectors
        __device__ void setDeviceRealValue(int r, int c, double val);       // Sets device real value for matrices
        __device__ void setDeviceImaginaryValue(int r, double val);           // Sets device imaginary value for vectors
        __device__ void setDeviceImaginaryValue(int r, int c, double val);    // Sets device imaginary value for matrices
        double getHostRealValue(int r);                                     // Gets host real value for vectors
        double getHostRealValue(int r, int c);                              // Gets host real value for matrices
        double getHostImaginaryValue(int r);                                  // Gets host imaginary value for vectors
        double getHostImaginaryValue(int r, int c);                           // Gets host imaginary value for matrices
        __device__ double getDeviceRealValue(int r);                        // Gets device real value for vectors
        __device__ double getDeviceRealValue(int r, int c);                 // Gets device real value for matrices
        __device__ double getDeviceImaginaryValue(int r);                     // Gets device imaginary value for vectors
        __device__ double getDeviceImaginaryValue(int r, int c);              // Gets device imaginary value for matrices

};


#endif