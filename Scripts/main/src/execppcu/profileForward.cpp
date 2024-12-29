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
    int iter;
    FILE *file;
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");
    char filename[256];
    if (deviceComputation) {
        sprintf(filename, "%s/Results/profiling/grid/forward_device.txt",dir);
        iter = 1e2;
    }
    else {
        iter = 10;
        int num_threads = omp_get_max_threads();
        printf("%d\n", num_threads);
        sprintf(filename, "%s/Results/profiling/grid/forward_host_%d.txt", dir, num_threads);
    }
    file = fopen(filename, "a");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points, deviceComputation, x, y, n, printOutput);
    
    bioscat.getNanostructure();

    // warm-up
    if (deviceComputation) {
        for (int i = 0; i < iter; i++) {
            bioscat.setupSegments();
        }
    }

    fprintf(file,"%d %d ", total_grid_points, n);
    double start_total = omp_get_wtime();
    double start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.setupSegments();
    }

    double end = omp_get_wtime();
    fprintf(file,"%e ",(end-start)/iter);

    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.reset();
    }
    end = omp_get_wtime();
    fprintf(file,"%e ",(end-start)/iter);

    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.prepareForward(beta, lambda);
    }
    end = omp_get_wtime();
    fprintf(file,"%e ",(end-start)/iter);

    for (int polarisation = 1; polarisation <= 2; polarisation++) {

        bioscat.polarisation = polarisation;
        bioscat.segments[0].polarisation = polarisation;
        
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.segments[0].computeFieldsForLinearSystem();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter); 
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.segments[0].setupRightHandSide();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter);
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.segments[0].setupSystemMatrix();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter);
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.segments[0].solveLinearSystem();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter);
       
        
        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.computeScatteredSubFields();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter);

        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.computeInteriorSubFields();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter);

        start = omp_get_wtime();
        for (int i = 0; i < iter; i++) {
            bioscat.computeIncidentSubFields();
        }
        end = omp_get_wtime();
        fprintf(file,"%e ",(end-start)/iter);

    }

    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.computeScatteredFields();
    }
    end = omp_get_wtime();
    fprintf(file,"%e ",(end-start)/iter);

     start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.computeInteriorFields();
    }
    end = omp_get_wtime();
    fprintf(file,"%e ",(end-start)/iter);


    start = omp_get_wtime();
    for (int i = 0; i < iter; i++) {
        bioscat.computeIncidentFields();
    }
    end = omp_get_wtime();
    fprintf(file,"%e ",(end-start)/iter);
    double end_total = omp_get_wtime();
    fprintf(file,"%e ",(end_total-start_total)/iter);
    bioscat.free();

    fprintf(file,"\n");

    fflush(file);
    fflush(stdout);
}



}