#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../../lib/BioScat.h"
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
extern "C" {
using namespace std;

void BioScat::forwardSolver(int polarisation) {

    this->polarisation = polarisation;

    double start, end, start_inner, end_inner;
    start = omp_get_wtime();
    //#pragma omp parallel for num_threads(num_segments)
    for (int i = 0; i < num_segments; i++) {
        segments[i].polarisation = polarisation;

        start_inner = omp_get_wtime();
        
        segments[i].computeFieldsForLinearSystem();
        segments[i].setupRightHandSide();
        segments[i].setupSystemMatrix();
        segments[i].solveLinearSystem();
      
        end_inner = omp_get_wtime();
        if (printOutput) printf("\nIt took %.4e seconds to solve the linear system for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }
    end = omp_get_wtime();
    if (printOutput) printf("\nIt took %.4e seconds to solve all the linear systems.\n\n",end - start);
}

}