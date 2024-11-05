#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
#include "../lib/BioScat.h"
using namespace std;


void generateArtificialData(double *x, double*y, int n, char * protein_structure, int num_segments, int total_grid_points, double * betas, double * lambdas, int num_betas, int num_lambdas) {
    RealMatrix reflectance = RealMatrix(n, num_betas, num_lambdas);
   
    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points);

    bioscat.getNanostructure();

    bioscat.getSegments();

    bioscat.setupObservationPoints(x, y, n);

    for (int i = 0; i < num_lambdas; i++) {
        
        double lambda = lambdas[i];
        bioscat.prepareForward(lambda);

        for (int polarisation = 1; polarisation <= 2; polarisation++) {
            bioscat.forwardSolver(polarisation);
            bioscat.computeScatteredSubFields();
            bioscat.computeReflectedSubFields();
            bioscat.computeIncidentSubFields();
        }

        for (int j = 0; j < num_betas; j++) {
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

    char filename[256];
    sprintf(filename, "../../../Data/artificial_data/temp/reflectance.txt");
    reflectance.dumpVector(filename);


}



}
