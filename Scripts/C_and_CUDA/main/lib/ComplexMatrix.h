#ifndef _COMPLEX_MATRIX_H
#define _COMPLEX_MATRIX_H

class ComplexMatrix {
    // Stores a complexMatrix rowmajor in a vector both on the host and the device
    private: // Private to avoid wrong indexing
        double * real_h;    // Real entries on host
        double * complex_h; // Complex entries on host
        double * real_d;    // Real entries on device
        double * complex_d; // Complex entries on device
    public:
        int    rows; // Number of rows
        int    cols; // Number of cols

        // Methods
        ComplexMatrix(int rows, int cols);                                  // Constructer, allocates arrays and initializes to zero
        void free();                                           // Frees arrays
        void toHost();                                                      // Sends data to host
        void toDevice();                                                    // Sends data to device
        void setHostRealValue(int r, double val);                           // Sets host real value for vectors
        void setHostRealValue(int r, int c, double val);                    // Sets host real value for matrices
        void setHostComplexValue(int r, double val);                        // Sets host complex value for vectors
        void setHostComplexValue(int r, int c, double val);                 // Sets host complex value for matrices
        __device__ void setDeviceRealValue(int r, double val);              // Sets device real value for vectors
        __device__ void setDeviceRealValue(int r, int c, double val);       // Sets device real value for matrices
        __device__ void setDeviceComplexValue(int r, double val);           // Sets device complex value for vectors
        __device__ void setDeviceComplexValue(int r, int c, double val);    // Sets device complex value for matrices
        double getHostRealValue(int r);                                     // Gets host real value for vectors
        double getHostRealValue(int r, int c);                              // Gets host real value for matrices
        double getHostComplexValue(int r);                                  // Gets host complex value for vectors
        double getHostComplexValue(int r, int c);                           // Gets host complex value for matrices
        __device__ double getDeviceRealValue(int r);                        // Gets device real value for vectors
        __device__ double getDeviceRealValue(int r, int c);                 // Gets device real value for matrices
        __device__ double getDeviceComplexValue(int r);                     // Gets device complex value for vectors
        __device__ double getDeviceComplexValue(int r, int c);              // Gets device complex value for matrices

};


#endif