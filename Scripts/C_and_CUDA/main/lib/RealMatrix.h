#ifndef _REAL_MATRIX_H
#define _REAL_MATRIX_H
extern "C" {
class RealMatrix {
    // Stores a real matrix rowmajor in a vector both on the host and the device
    public:
        // Variables
        double * val_h; // Entries on host
        double * val_d; // Entries on device
        int    rows; // Number of rows
        int    cols; // Number of cols

        // Methods
        RealMatrix() {
            rows = 0;
            cols = 0;
            val_h = NULL;
            val_d = NULL;
        };  // Constructer
        RealMatrix(int rows, int cols);  // Constructer, allocates arrays and initializes to zero
        void freeRealMatrix();                   // Destructer, frees arrays
        void toHost();                             // Sends data to host
        void toDevice();                           // Sends data to device

};
}

#endif