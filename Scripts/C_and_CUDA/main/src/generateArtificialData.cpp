#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
#include "../lib/BioScat.h"
using namespace std;

void showProgressBar(float progress) {
    // From chatgpt
    int barWidth = 100; // Width of the progress bar in characters

    std::cout << " ";
    int pos = barWidth * progress; // Position of the current progress in the bar
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "|";
        //else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << " " << int(progress * 100.0) << " %\r"; // Display percentage
    std::cout.flush();
}


void generateArtificialData(double *x, double*y, int n, char * protein_structure, int num_segments, int total_grid_points, double * betas, double * lambdas, int num_betas, int num_lambdas) {
    
    RealMatrix reflectance = RealMatrix(n, num_betas, num_lambdas);
   

    double start = omp_get_wtime();
    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points);
    bioscat.printOutput = false;
    
    bioscat.setupObservationPoints(x, y, n);

    bioscat.getNanostructure();

    bioscat.getSegments();
    
    
    for (int i = 0; i < num_lambdas; i++) {
        bioscat.reset();
        double lambda = lambdas[i];
        bioscat.prepareForward(lambda);

        for (int polarisation = 1; polarisation <= 2; polarisation++) {
            bioscat.forwardSolver(polarisation);
            bioscat.computeScatteredSubFields();
            bioscat.computeReflectedSubFields();
            bioscat.computeIncidentSubFields();
        }

        for (int j = 0; j < num_betas; j++) {
            showProgressBar((i*num_betas + j) / ((double) num_lambdas*num_betas));
            double beta = betas[j];

            bioscat.computeScatteredFields(beta);

            bioscat.computeReflectedFields(beta);

            bioscat.computeIncidentFields(beta);

            bioscat.computeReflectance();

            for (int k = 0; k < n; k++) {
                reflectance.setHostValue(k,i,j,bioscat.reflectance.getHostValue(k));
            }
        }
    }
    showProgressBar(1.0);
    double end = omp_get_wtime();
    printf("\nIt took %.4f to compute the reflectance!\n",end-start);

    char filename[256];
    sprintf(filename, "../../../Data/artificial_data/temp/reflectance.txt");
    reflectance.dumpVector(filename);

    bioscat.free();


}



}
