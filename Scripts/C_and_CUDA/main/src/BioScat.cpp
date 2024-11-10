
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include <cblas.h>
#include <math.h>
extern "C" {
#include "../lib/BioScat.h"
#include "../lib/Segment.h"
#include "../lib/RealMatrix.h"
#include "../lib/combinePolarisation.h"
using namespace std;

BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points) {

    this->protein_structure = protein_structure;
    this->num_segments = num_segments;
    this->total_grid_points = total_grid_points;

}

void BioScat::free() {

    for (int i = 0; i < num_segments; i++) {
        segments[i].free();
    }

    delete[] segments;
    reflectance.free();
    x_obs.free();
    y_obs.free();

    E_scat.free();
    H_scat.free();
    E_inc.free();
    H_inc.free();
    E_ref.free();
    H_ref.free();
    for (int i = 0; i < 2; i++) {
        E_scat_pol[i].free();
        H_scat_pol[i].free();
        E_inc_pol[i].free();
        H_inc_pol[i].free();
        E_ref_pol[i].free();
        H_ref_pol[i].free();
    }
    //cudaDeviceReset();
    printf("DESTRUCTED!\n");

}

void BioScat::getSegments() {

    if (deviceComputation) {
        if (printOutput) printf("We are computing on the device!\n");
    }
    else {
        if (printOutput) printf("We are computing on the host!\n");
    }

    double start = omp_get_wtime();

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].deviceComputation = deviceComputation;
        segments[i].current_segment = i;
        segments[i].n_obs = n_obs;
        segments[i].setup(this->nanostructure, total_grid_points, num_segments);
    }


    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %e seconds to setup the segments!\n\n",end-start);

}

void BioScat::getSegments(Nanostructure nanostructure) {

    double start = omp_get_wtime();

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        //segments[i] = Segment();
        segments[i].deviceComputation = deviceComputation;
        segments[i].current_segment = i;
        segments[i].n_obs = n_obs;
        segments[i].setup(nanostructure, total_grid_points, num_segments);
    }

    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %e seconds to setup the segments!\n\n",end-start);
}

void BioScat::prepareForward(double lambda) {
    for (int i = 0; i < num_segments; i++) {
        segments[i].newWavelength(lambda);
    }
}

void BioScat::prepareForward(double beta, double lambda) {
    this->beta = beta;
    for (int i = 0; i < num_segments; i++) {
        segments[i].newWavelength(lambda);
    }
}

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
        //segments[i].freeScatteredFields();
        //segments[i].freeInteriorFields();
        //segments[i].freeIncidentFields();
        //segments[i].freeReflectedFields();
        segments[i].solveLinearSystem();
        end_inner = omp_get_wtime();
        if (printOutput) printf("\nIt took %.4e seconds to solve the linear system for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }
    end = omp_get_wtime();
    if (printOutput) printf("\nIt took %.4e seconds to solve all the linear systems.\n\n",end - start);

}

void BioScat::setupObservationPoints(double *x, double*y, int n) {
    n_obs = n;
    x_obs = RealMatrix(n);
    y_obs = RealMatrix(n);
    for (int i = 0; i < n; i++) x_obs.setHostValue(i, x[i]);
    for (int i = 0; i < n; i++) y_obs.setHostValue(i, y[i]);
    x_obs.toDevice();
    y_obs.toDevice();
    /*x_obs.rows = n;
    y_obs.rows = n;
    x_obs.cols = 1;
    y_obs.cols = 1;
    x_obs.setHostPointer(x);
    y_obs.setHostPointer(y);*/

    // Allocate fields
    E_scat = Field(n);
    H_scat = Field(n);
    E_inc  = Field(n);
    H_inc  = Field(n);
    E_ref  = Field(n);
    H_ref  = Field(n);
    for (int i = 0; i < 2; i++) {
        E_scat_pol[i] = Field(n);
        H_scat_pol[i] = Field(n);
        E_inc_pol[i]  = Field(n);
        H_inc_pol[i]  = Field(n);
        E_ref_pol[i]  = Field(n);
        H_ref_pol[i]  = Field(n);
    }

    // Allocate reflectance array
    reflectance = RealMatrix(n);
}

void BioScat::computeScatteredFields() {

    computeScatteredFields(beta);
}


void BioScat::computeScatteredFields(double beta) {
    combinePolarisation(E_scat_pol, E_scat, beta);
    combinePolarisation(H_scat_pol, H_scat, beta);
    
}

void BioScat::computeIncidentFields() {
    
    computeIncidentFields(beta);
}


void BioScat::computeIncidentFields(double beta) {

    combinePolarisation(E_inc_pol, E_inc, beta);
    combinePolarisation(H_inc_pol, H_inc, beta);

    E_inc.toHost();
}

void BioScat::computeReflectedFields() {
    computeReflectedFields(beta);
}


void BioScat::computeReflectedFields(double beta) {

    combinePolarisation(E_ref_pol, E_ref, beta);
    combinePolarisation(H_ref_pol, H_ref, beta);

    E_ref.toHost();

    
}

void BioScat::dumpFields() {

    char filename[256];

    // Save scattered electric fields
    sprintf(filename, "../../../Results/forward/Ex_scat.txt");
    E_scat.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ey_scat.txt");
    E_scat.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ez_scat.txt");
    E_scat.z.dumpResult(filename);

    // Save scattered magnetic fields
    sprintf(filename, "../../../Results/forward/Hx_scat.txt");
    H_scat.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hy_scat.txt");
    H_scat.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hz_scat.txt");
    H_scat.z.dumpResult(filename);

    // Save incident electric fields
    sprintf(filename,"../../../Results/forward/Ex_inc.txt");
    E_inc.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ey_inc.txt");
    E_inc.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ez_inc.txt");
    E_inc.z.dumpResult(filename);

    // Save incident magnetic fields
    sprintf(filename,"../../../Results/forward/Hx_inc.txt");
    H_inc.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hy_inc.txt");
    H_inc.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hz_inc.txt");
    H_inc.z.dumpResult(filename);

     // Save reflected electric fields
    sprintf(filename,"../../../Results/forward/Ex_ref.txt");
    E_ref.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ey_ref.txt");
    E_ref.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Ez_ref.txt");
    E_ref.z.dumpResult(filename);

    // Save reflected magnetic fields
    sprintf(filename,"../../../Results/forward/Hx_ref.txt");
    H_ref.x.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hy_ref.txt");
    H_ref.y.dumpResult(filename);
    sprintf(filename,"../../../Results/forward/Hz_ref.txt");
    H_ref.z.dumpResult(filename);
}

double getSquaredNorm(double x_real, double x_imag, double y_real, double y_imag, double z_real, double z_imag) {
    double x2 = x_real*x_real + x_imag*x_imag;
    double y2 = y_real*y_real + y_imag*y_imag;
    double z2 = z_real*z_real + z_imag*z_imag;
    return x2 + y2 + z2;
}

void BioScat::computeReflectance() {
    //#pragma omp parallel for num_threads(x_obs.rows)
    for (int i = 0; i < x_obs.rows; i++) {
        double numerator = getSquaredNorm(E_inc.x.getHostRealValue(i), E_inc.x.getHostImagValue(i),
                                          E_inc.y.getHostRealValue(i), E_inc.y.getHostImagValue(i),
                                          E_inc.z.getHostRealValue(i), E_inc.z.getHostImagValue(i));
        double denominator = getSquaredNorm(E_ref.x.getHostRealValue(i) + E_scat.x.getHostRealValue(i), E_ref.x.getHostImagValue(i) + E_scat.x.getHostImagValue(i),
                                            E_ref.y.getHostRealValue(i) + E_scat.y.getHostRealValue(i), E_ref.y.getHostImagValue(i) + E_scat.y.getHostImagValue(i),
                                            E_ref.z.getHostRealValue(i) + E_scat.z.getHostRealValue(i), E_ref.z.getHostImagValue(i) + E_scat.z.getHostImagValue(i));
        reflectance.setHostValue(i, numerator/denominator);
    }
}

void BioScat::reset() {
    if (true) {
        // Initialize all arrays to zero on the host
        for (int i = 0; i < 2; i++) {
            E_scat_pol[i].setDeviceZero();
            H_scat_pol[i].setDeviceZero();
            E_inc_pol[i].setDeviceZero();
            H_inc_pol[i].setDeviceZero();
            E_ref_pol[i].setDeviceZero();
            H_ref_pol[i].setDeviceZero();
        }
    }
    else {
        // Initialize all arrays to zero on the host
        for (int i = 0; i < 2; i++) {
            E_scat_pol[i].setHostZero();
            H_scat_pol[i].setHostZero();
            E_inc_pol[i].setHostZero();
            H_inc_pol[i].setHostZero();
            E_ref_pol[i].setHostZero();
            H_ref_pol[i].setHostZero();
        }
    }
}

}