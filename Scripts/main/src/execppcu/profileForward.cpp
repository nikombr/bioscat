#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "../../lib/utils/ComplexMatrix.h"
#include "../../lib/BioScat.h"
#include <omp.h>
extern "C" {
using namespace std;


void profileForward(double *x, double*y, int n ,char * protein_structure, int num_segments, int total_grid_points, double beta, double lambda, int deviceComputation_int) {
    bool deviceComputation = deviceComputation_int == 1 ? true : false;
    bool printOutput = false;

    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points, deviceComputation, x, y, n, printOutput);
    
    bioscat.getNanostructure();

    int iter = 1e2;
    double start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.setupSegments();
    }

    double end = omp_get_wtime();
    printf("%e ",(end-start)/iter);

    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.reset();
    }
    end = omp_get_wtime();
    printf("%e ",(end-start)/iter);

    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.prepareForward(beta, lambda);
    }
    end = omp_get_wtime();
    printf("%e ",(end-start)/iter);


    for (int polarisation = 1; polarisation <= 2; polarisation++) {
       
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.forwardSolver(polarisation);
        }
        end = omp_get_wtime();
        printf("%e ",(end-start)/iter); 
        
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.computeScatteredSubFields();
        }
        end = omp_get_wtime();
        printf("%e ",(end-start)/iter);

        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.computeInteriorSubFields();
        }
        end = omp_get_wtime();
        printf("%e ",(end-start)/iter);

        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.computeIncidentSubFields();
        }
        end = omp_get_wtime();
        printf("%e ",(end-start)/iter);

    }

    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.computeScatteredFields();
    }
    end = omp_get_wtime();
    printf("%e ",(end-start)/iter);

     start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.computeInteriorFields();
    }
    end = omp_get_wtime();
    printf("%e ",(end-start)/iter);


    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.computeIncidentFields();
    }
    end = omp_get_wtime();
    printf("%e ",(end-start)/iter);
    
    bioscat.free();

    printf("\n");

}



}