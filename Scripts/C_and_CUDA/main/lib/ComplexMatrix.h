#ifndef _COMPLEX_MATRIX_H
#define _COMPLEX_MATRIX_H

class ComplexMatrix {
    // Stores a complexMatrix rowmajor in a vector both on the host and the device
    public:
        // Variables
        double * real_h; // Real entries on host
        double * complex_h; // Complex entries on host
        double * real_d; // Real entries on device
        double * complex_d; // Complex entries on device
        int    rows; // Number of rows
        int    cols; // Number of cols

        // Methods
        ComplexMatrix(int rows, int cols);  // Constructer, allocates arrays and initializes to zero
        ~ComplexMatrix();                   // Destructer, frees arrays
        void toHost();                             // Sends data to host
        void toDevice();                           // Sends data to device

};


#endif