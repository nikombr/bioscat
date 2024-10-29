
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include <cblas.h>
extern "C" {
#include "../lib/BioScat.h"
#include "../lib/Segment.h"
#include "../lib/RealMatrix.h"
using namespace std;

BioScat::BioScat(char* protein_structure, int num_segments) {

    this->protein_structure = protein_structure;
    this->num_segments = num_segments;

}

BioScat::~BioScat() {

    for (int i = 0; i < num_segments; i++) {
        segments[i].free();
    }

    delete[] segments;

    printf("DESTRUCTED!\n");

}

void BioScat::getSegments() {

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].setup(this->nanostructure, i, total_grid_points, num_segments);
    }

}

void BioScat::getSegments(Nanostructure nanostructure) {

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].setup(nanostructure, i, total_grid_points, num_segments);
    }

}

void BioScat::forwardSolver(int scenario) {

    double start, end;

    for (int i = 0; i < num_segments; i++) {
        segments[i].scenario = scenario;

        start = omp_get_wtime();
        segments[i].computeFieldsForLinearSystem();
        segments[i].setupRightHandSide();
        segments[i].setupSystemMatrix();
        segments[i].solveLinearSystem();
        end = omp_get_wtime();
        printf("\nIt took %.4e seconds to solve the linear system for segment %d.\n\n",end-start,i+1);
        
    }

}

void BioScat::setupObservationPoints(double *x, double*y, int n) {
    x_obs = RealMatrix();
    y_obs = RealMatrix();
    x_obs.rows = n;
    y_obs.rows = n;
    x_obs.cols = 1;
    y_obs.cols = 1;
    x_obs.setHostPointer(x);
    y_obs.setHostPointer(y);
}

void BioScat::computeScatteredFields() {

    int n = x_obs.rows;
    double val;

    // Allocate fields (we may need to initialize to zero)
    E_scat = Field(n, true, true, true);
    H_scat = Field(n, true, true, true);

    double start = omp_get_wtime();
    for (int i = 0; i < num_segments; i++) {

        
        segments[i].computeScatteredFieldMatrices(x_obs, y_obs, false);


        if (segments[i].scenario == 1) {
            for (int k = 0; k < n; k++) {
                val = E_scat.z.getHostRealValue(k);
                for (int j = 0; j < segments[i].n_int; j++) {
                    val += segments[i].E_scat_matrix.z.getHostRealValue(k,j) * segments[i].C.getHostRealValue(j);
                    val -= segments[i].E_scat_matrix.z.getHostImagValue(k,j) * segments[i].C.getHostImagValue(j);
                }
                E_scat.z.setHostRealValue(k, val);
            }
            for (int k = 0; k < n; k++) {
                val = E_scat.z.getHostImagValue(k);
                for (int j = 0; j < segments[i].n_int; j++) {
                    val += segments[i].E_scat_matrix.z.getHostRealValue(k,j) * segments[i].C.getHostImagValue(j);
                    val += segments[i].E_scat_matrix.z.getHostImagValue(k,j) * segments[i].C.getHostRealValue(j);
                }
                E_scat.z.setHostImagValue(k, val);
            }

            for (int k = 0; k < n; k++) {
                val = H_scat.x.getHostRealValue(k);
                for (int j = 0; j < segments[i].n_int; j++) {
                    val += segments[i].H_scat_matrix.x.getHostRealValue(k,j) * segments[i].C.getHostRealValue(j);
                    val -= segments[i].H_scat_matrix.x.getHostImagValue(k,j) * segments[i].C.getHostImagValue(j);
                }
                H_scat.x.setHostRealValue(k, val);
            }
            for (int k = 0; k < n; k++) {
                val = H_scat.x.getHostImagValue(k);
                for (int j = 0; j < segments[i].n_int; j++) {
                    val += segments[i].H_scat_matrix.x.getHostRealValue(k,j) * segments[i].C.getHostImagValue(j);
                    val += segments[i].H_scat_matrix.x.getHostImagValue(k,j) * segments[i].C.getHostRealValue(j);
                }
                H_scat.x.setHostImagValue(k, val);
            }

            for (int k = 0; k < n; k++) {
                val = H_scat.y.getHostRealValue(k);
                for (int j = 0; j < segments[i].n_int; j++) {
                    val += segments[i].H_scat_matrix.y.getHostRealValue(k,j) * segments[i].C.getHostRealValue(j);
                    val -= segments[i].H_scat_matrix.y.getHostImagValue(k,j) * segments[i].C.getHostImagValue(j);
                }
                H_scat.y.setHostRealValue(k, val);
            }
            for (int k = 0; k < n; k++) {
                val = H_scat.y.getHostImagValue(k);
                for (int j = 0; j < segments[i].n_int; j++) {
                    val += segments[i].H_scat_matrix.y.getHostRealValue(k,j) * segments[i].C.getHostImagValue(j);
                    val += segments[i].H_scat_matrix.y.getHostImagValue(k,j) * segments[i].C.getHostRealValue(j);
                }
                H_scat.y.setHostImagValue(k, val);
            }

        }
        printf("hmm");
        
        //segments[i].computeInteriorFieldMatrices(x_obs, y_obs);
        
    }

    double end = omp_get_wtime();
    printf("\nIt took %.4e seconds to compute the scattered fields in the observation points.\n\n",end-start);

    FILE *file;
    char * filename = "../../../Results/forward/Ez_scat.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < n; r++) {
        fprintf(file, "%e\t%e\n", E_scat.z.getHostRealValue(r), E_scat.z.getHostImagValue(r));
    }
    fclose(file);

    filename = "../../../Results/forward/Hx_scat.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < n; r++) {
        fprintf(file, "%e\t%e\n", H_scat.x.getHostRealValue(r), H_scat.x.getHostImagValue(r));
    }
    fclose(file);
    filename = "../../../Results/forward/Hy_scat.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < n; r++) {
        fprintf(file, "%e\t%e\n", H_scat.y.getHostRealValue(r), H_scat.y.getHostImagValue(r));
    }
    fclose(file);
}


void BioScat::computeIncidentFields() {

    int n = x_obs.rows;
    double val;

    // Allocate fields (we may need to initialize to zero)
    E_inc = Field(n, true, true, true);
    H_inc = Field(n, true, true, true);

    double start = omp_get_wtime();

    segments[0].computeIncidentFieldVectors(y_obs);

    if (segments[0].scenario == 1) {
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.z.getHostRealValue(k);
            E_inc.z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.z.getHostImagValue(k);
            E_inc.z.setHostImagValue(k, val);
        }
    }


    for (int i = 1; i < num_segments; i++) {

        segments[i].computeIncidentFieldVectors(y_obs);

        if (segments[i].scenario == 1) {
            for (int k = 0; k < n; k++) {
                E_inc.z.setHostRealValue(k, segments[i].E_inc_vector.z.getHostRealValue(k));
            }
            for (int k = 0; k < n; k++) {
                E_inc.z.setHostImagValue(k, segments[i].E_inc_vector.z.getHostImagValue(k));
            }
            

        }
        
        //segments[i].computeInteriorFieldMatrices(x_obs, y_obs);
        
    }

    double end = omp_get_wtime();
    printf("\nIt took %.4e seconds to compute the incident fields in the observation points.\n\n",end-start);

    FILE *file;
    char * filename = "../../../Results/forward/Ez_inc.txt";
    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int r = 0; r < n; r++) {
        fprintf(file, "%e\t%e\n", E_inc.z.getHostRealValue(r), E_inc.z.getHostImagValue(r));
    }
    fclose(file);
}

void BioScat::computeReflectedFields() {
    int j = 0;
}

}