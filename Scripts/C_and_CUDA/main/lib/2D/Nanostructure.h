#ifndef _NANOSTRUCTURE_H
#define _NANOSTRUCTURE_H
extern "C" {
#include "../RealMatrix.h"

struct Nanostructure {
    RealMatrix f; // Function value in x
    RealMatrix x; // x-values

    // Constructor 
    Nanostructure() {
        f = RealMatrix();
        x = RealMatrix();
        std::cout << "Empty segment constructor." << std::endl;
    } 
    Nanostructure(int n) {
        f = RealMatrix(n, 1);
        x = RealMatrix(n, 1);
        std::cout << "Non-empty segment constructor." << std::endl;
    }
};
}

#endif