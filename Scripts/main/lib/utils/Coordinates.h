#ifndef _COORDINATES_H
#define _COORDINATES_H
#include "utils/RealMatrix.h"
extern "C" {

struct Coordinates {
    RealMatrix x; // Function value in x
    RealMatrix y; // x-values


    // Constructor 
    Coordinates() {
        x = RealMatrix();
        y = RealMatrix();
        //std::cout << "Empty segment constructor." << std::endl;
    } 
    Coordinates(int n) {
        x = RealMatrix(n);
        y = RealMatrix(n);
        //std::cout << "Non-empty segment constructor." << std::endl;
    }

    Coordinates(int n, bool host, bool device) {
        x = RealMatrix(n, host, device);
        y = RealMatrix(n, host, device);
        //std::cout << "Non-empty segment constructor." << std::endl;
    }

    void free() {
        x.free();
        y.free();
    }

    void allocateHost() {
        x.allocateHost();
        y.allocateHost();
    }

    void toHost() {
        x.toHost();
        y.toHost();
    }
};


}

#endif