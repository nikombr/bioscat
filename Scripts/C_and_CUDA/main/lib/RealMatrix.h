#ifndef _REAL_MATRIX_H
#define _REAL_MATRIX_H
extern "C" {
class RealMatrix {
    // Stores a real matrix rowmajor in a vector both on the host and the device
    private: // Private to avoid wrong indexing
        double * val_h; // Entries on host
        double * val_d; // Entries on device
    public:
        int rows; // Number of rows
        int cols; // Number of cols

        // Methods
        RealMatrix() {
            rows = 0;
            cols = 0;
            val_h = NULL;
            val_d = NULL;
        };  // Constructer
        RealMatrix(int rows);                                       // Constructer, allocates vectors and initializes to zero
        RealMatrix(int rows, int cols);                             // Constructer, matrices arrays and initializes to zero
        void free();                                                // Frees arrays
        void toHost();                                              // Sends data from device to host
        void toDevice();                                            // Sends data from host to device
        void setHostValue(int r, double val);                       // Sets host value for vectors
        void setHostValue(int r, int c, double val);                // Sets host value for matrices
        __device__ void setDeviceValue(int r, double val) {         // Sets device value for vectors
            val_d[r] = val;
        }          
        __device__ void setDeviceValue(int r, int c, double val) {  // Sets device value for matrices
            val_d[r*cols + c] = val;
        }   
        double getHostValue(int r);                                 // Gets host value for vectors
        double getHostValue(int r, int c);                          // Gets host value for matrices
        __device__ double getDeviceValue(int r) {                   // Gets device value for vectors
            return val_d[r];
        }                    
        __device__ double getDeviceValue(int r, int c) {            // Gets device value for matrices
            return val_d[r*cols + c];
        }            
        double * getHostPointer() {
            return val_h;
        }
        double * getDevicePointer() {
            return val_d;
        }
        void setHostPointer(double * val) {
            val_h = val;
        }
   

};
}

#endif