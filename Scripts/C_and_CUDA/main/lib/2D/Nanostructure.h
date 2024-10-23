#ifndef _NANOSTRUCTURE_H
#define _NANOSTRUCTURE_H

#include "../RealMatrix.h"

struct Nanostructure {
    RealMatrix f; // Function value in x
    RealMatrix x; // x-values

    // Constructor 
    Nanostructure(int n) : f(n, 1), x(n, 1) {
        
    }
};

#endif