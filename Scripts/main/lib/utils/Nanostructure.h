#ifndef _NANOSTRUCTURE_H
#define _NANOSTRUCTURE_H
#include "utils/RealMatrix.h"
extern "C" {

struct Nanostructure {
    RealMatrix f; // Function value in x
    RealMatrix x; // x-values

    // Constructor 
    Nanostructure() {
        f = RealMatrix();
        x = RealMatrix();
        //std::cout << "Empty segment constructor." << std::endl;
    } 
    Nanostructure(int n) {
        f = RealMatrix(n, 1);
        x = RealMatrix(n, 1);
        //std::cout << "Non-empty segment constructor." << std::endl;
    }

    void free() {
        f.free();
        x.free();
    }
};


}

#endif