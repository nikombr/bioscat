
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

    if (deviceComputation) {
        printf("We are computing on the device!\n");
    }
    else {
        printf("We are computing on the host!\n");
    }

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].deviceComputation = deviceComputation;
        segments[i].current_segment = i;
        segments[i].setup(this->nanostructure, total_grid_points, num_segments);
    }

    printf("Segments are ready!\n");

}

void BioScat::getSegments(Nanostructure nanostructure) {

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i] = Segment();
        segments[i].deviceComputation = deviceComputation;
        segments[i].current_segment = i;
        segments[i].setup(nanostructure, total_grid_points, num_segments);
    }

}

void BioScat::forwardSolver(int polarisation) {

    this->polarisation = polarisation;

    double start, end, start_inner, end_inner;
    start = omp_get_wtime();
    for (int i = 0; i < num_segments; i++) {
        segments[i].polarisation = polarisation;

        start_inner = omp_get_wtime();
        segments[i].computeFieldsForLinearSystem();
        segments[i].setupRightHandSide();
        segments[i].setupSystemMatrix();
        segments[i].solveLinearSystem();
        end_inner = omp_get_wtime();
        printf("\nIt took %.4e seconds to solve the linear system for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }
    end = omp_get_wtime();
    printf("\nIt took %.4e seconds to solve all the linear systems.\n\n",end - start);

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
    

    // Allocate fields
    E_scat[polarisation - 1] = Field(n, true, true, true);
    H_scat[polarisation - 1] = Field(n, true, true, true);
    
    double start = omp_get_wtime();
    for (int i = 0; i < num_segments; i++) {
        
        double start_inner = omp_get_wtime();
        
        segments[i].computeScatteredFieldMatrices(x_obs, y_obs, false);

        double end_inner = omp_get_wtime();
        printf("\nIt took %.4e seconds to compute the scattered field matrices for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }
    
    if (polarisation == 1) {
        #pragma omp parallel for
        for (int k = 0; k < n; k++) {
            double val1, val2, val3, C_real, C_imag;
            val1 = 0.0;
            val2 = 0.0;
            val3 = 0.0;
            for (int i = 0; i < num_segments; i++) {
                
                for (int j = 0; j < segments[i].n_int; j++) {
                    C_real = segments[i].C.getHostRealValue(j);
                    C_imag = segments[i].C.getHostImagValue(j);
                    val1 += segments[i].E_scat_matrix.z.getHostRealValue(k,j) * C_real;
                    val1 -= segments[i].E_scat_matrix.z.getHostImagValue(k,j) * C_imag;
                    val2 += segments[i].H_scat_matrix.x.getHostRealValue(k,j) * C_real;
                    val2 -= segments[i].H_scat_matrix.x.getHostImagValue(k,j) * C_imag;
                    val3 += segments[i].H_scat_matrix.y.getHostRealValue(k,j) * C_real;
                    val3 -= segments[i].H_scat_matrix.y.getHostImagValue(k,j) * C_imag;
            
                }
                
            }

            E_scat[0].z.setHostRealValue(k, val1);
            H_scat[0].x.setHostRealValue(k, val2);
            H_scat[0].y.setHostRealValue(k, val3);
        }
    
        #pragma omp parallel for
        for (int k = 0; k < n; k++) {
            double val1, val2, val3, C_real, C_imag;
            val1 = 0.0;
            val2 = 0.0;
            val3 = 0.0;
            for (int i = 0; i < num_segments; i++) {
                //if (i == 1) {
                for (int j = 0; j < segments[i].n_int; j++) {
                    C_real = segments[i].C.getHostRealValue(j);
                    C_imag = segments[i].C.getHostImagValue(j);
                    val1 += segments[i].E_scat_matrix.z.getHostRealValue(k,j) * C_imag;
                    val1 += segments[i].E_scat_matrix.z.getHostImagValue(k,j) * C_real;
                    val2 += segments[i].H_scat_matrix.x.getHostRealValue(k,j) * C_imag;
                    val2 += segments[i].H_scat_matrix.x.getHostImagValue(k,j) * C_real;
                    val3 += segments[i].H_scat_matrix.y.getHostRealValue(k,j) * C_imag;
                    val3 += segments[i].H_scat_matrix.y.getHostImagValue(k,j) * C_real;
                }
                //}
            }
            E_scat[0].z.setHostImagValue(k, val1);
            H_scat[0].x.setHostImagValue(k, val2);
            H_scat[0].y.setHostImagValue(k, val3);
            
        }

        

           
    }
    double end = omp_get_wtime();
    printf("\nIt took %.4e seconds to compute the combined scattered fields in the observation points.\n\n",end-start);

}


void BioScat::computeIncidentFields() {

    int n = x_obs.rows;
    double val;

    // Allocate fields (we may need to initialize to zero)
    E_inc[polarisation - 1] = Field(n, true, true, true);
    H_inc[polarisation - 1] = Field(n, true, true, true);

    double start = omp_get_wtime();

    segments[0].computeIncidentFieldVectors(y_obs);

    if (segments[0].polarisation == 1) {
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.z.getHostRealValue(k);
            E_inc[0].z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.z.getHostImagValue(k);
            E_inc[0].z.setHostImagValue(k, val);
        }

        for (int k = 0; k < n; k++) {
            val = segments[0].H_inc_vector.x.getHostRealValue(k);
            H_inc[0].x.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].H_inc_vector.x.getHostImagValue(k);
            H_inc[0].x.setHostImagValue(k, val);
        }
    }
}

void BioScat::computeReflectedFields() {

    int n = x_obs.rows;
    double val;

    // Allocate fields (we may need to initialize to zero)
    E_ref[polarisation - 1] = Field(n, true, true, true);
    H_ref[polarisation - 1] = Field(n, true, true, true);

    double start = omp_get_wtime();

    segments[0].computeIncidentFieldVectors(y_obs);

    if (segments[0].polarisation == 1) {
        for (int k = 0; k < n; k++) {
            val = segments[0].E_ref_vector.z.getHostRealValue(k);
            E_ref[0].z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_ref_vector.z.getHostImagValue(k);
            E_ref[0].z.setHostImagValue(k, val);
        }

        for (int k = 0; k < n; k++) {
            val = segments[0].H_ref_vector.x.getHostRealValue(k);
            H_ref[0].x.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].H_ref_vector.x.getHostImagValue(k);
            H_ref[0].x.setHostImagValue(k, val);
        }
    }
}

void BioScat::dumpFields() {

    char filename[256];

    if (polarisation == 1) {

        // Save non-zero scattered fields
        sprintf(filename,"../../../Results/forward/Ez_scat_polarisation_%d.txt",polarisation);
        E_scat[0].z.dumpResult(filename);
        sprintf(filename, "../../../Results/forward/Hx_scat_polarisation_%d.txt",polarisation);
        H_scat[0].x.dumpResult(filename);
        sprintf(filename,"../../../Results/forward/Hy_scat_polarisation_%d.txt",polarisation);
        H_scat[0].y.dumpResult(filename);

        // Save non-zero interior fields
        /*filename = "../../../Results/forward/Ez_int.txt";
        E_int.z.dumpResult(filename.c_str());
        filename = "../../../Results/forward/Hx_int.txt";
        H_int.x.dumpResult(filename.c_str());
        filename = "../../../Results/forward/Hy_int.txt";
        H_int.y.dumpResult(filename.c_str());*/

        // Save non-zero incident fields
        sprintf(filename,"../../../Results/forward/Ez_inc_polarisation_%d.txt",polarisation);
        E_inc[0].z.dumpResult(filename);
        sprintf(filename,"../../../Results/forward/Hx_inc_polarisation_%d.txt",polarisation);
        H_inc[0].x.dumpResult(filename);

        // Save non-zero reflected fields
        sprintf(filename,"../../../Results/forward/Ez_ref_polarisation_%d.txt",polarisation);
        E_ref[0].z.dumpResult(filename);
        sprintf(filename,"../../../Results/forward/Hx_ref_polarisation_%d.txt",polarisation);
        H_ref[0].x.dumpResult(filename);

    }
    else if (polarisation == 2) {

        // Save non-zero scattered fields
        sprintf(filename,"../../../Results/forward/Hz_scat_polarisation_%d.txt",polarisation);
        H_scat[1].z.dumpResult(filename);
        sprintf(filename, "../../../Results/forward/Ex_scat_polarisation_%d.txt",polarisation);
        E_scat[1].x.dumpResult(filename);
        sprintf(filename,"../../../Results/forward/Ey_scat_polarisation_%d.txt",polarisation);
        E_scat[1].y.dumpResult(filename);

        // Save non-zero incident fields
        sprintf(filename,"../../../Results/forward/Hz_inc_polarisation_%d.txt",polarisation);
        H_inc[1].z.dumpResult(filename);
        sprintf(filename,"../../../Results/forward/Ex_inc_polarisation_%d.txt",polarisation);
        E_inc[1].x.dumpResult(filename);

        // Save non-zero reflected fields
        sprintf(filename,"../../../Results/forward/Hz_ref_polarisation_%d.txt",polarisation);
        H_ref[1].z.dumpResult(filename);
        sprintf(filename,"../../../Results/forward/Ex_ref_polarisation_%d.txt",polarisation);
        E_ref[1].x.dumpResult(filename);

    }

    
}

}