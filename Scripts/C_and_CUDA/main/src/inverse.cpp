#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include <math.h>
extern "C" {
#include "../lib/ComplexMatrix.h"
#include "../lib/BioScat.h"
using namespace std;

int countNumbersInFile(char* filename) {
    double num;
    FILE * file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file: %s\n",filename);
        return 0;
    }

    int count = 0;

    while (fscanf(file, "%e", &num) == 1) {
        count++;
    }

    // Close the file
    fclose(file);

    return count;
}

void getNumPoints(int * num_obs_points, int * num_lambdas, int * num_betas, char * protein_structure, int total_grid_points) {
    // Get number of data points
    char filename[256];

    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/x_obs.txt",protein_structure,total_grid_points);

    * num_obs_points = countNumbersInFile(filename);

    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/lambdas.txt",protein_structure,total_grid_points);

    * num_lambdas = countNumbersInFile(filename);
   
    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/betas.txt",protein_structure,total_grid_points);

    * num_betas = countNumbersInFile(filename);
    
}

void loadData(RealMatrix trueReflectance, RealMatrix lambdas, RealMatrix betas, RealMatrix x_obs, RealMatrix y_obs, char * protein_structure, int total_grid_points) {
    char filename[256];
    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/reflectance.txt",protein_structure,total_grid_points);
    trueReflectance.loadVector(filename);
    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/lambdas.txt",protein_structure,total_grid_points);
    lambdas.loadVector(filename);
    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/betas.txt",protein_structure,total_grid_points);
    betas.loadVector(filename);
    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/x_obs.txt",protein_structure,total_grid_points);
    x_obs.loadVector(filename);
    sprintf(filename, "../../../Data/artificial_data/%s/num_segments_1_total_grid_points_%d/y_obs.txt",protein_structure,total_grid_points);
    y_obs.loadVector(filename);
}

void initCoordInNanostructure(Nanostructure nanostructure, char * protein_structure, int total_grid_points) {
    double val;
    char filename[256];
    sprintf(filename, "../../../Data/nanostructures/%s_2D_x_%d.txt", protein_structure, total_grid_points);
    //printf("filename = %s\n",filename);

    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        printf("File: %s\n",filename);
        return;
    }
    for (int i = 0; i < total_grid_points; i++) {
        fscanf(file, "%lf,", &val);  // Reading each value into the array
        nanostructure.x.setHostValue(i, val);
    }
    fclose(file);
}

double computeLogLikelihood(RealMatrix trueReflectance, RealMatrix reflectance) {
    double beta = 1e6;
    double L = 0;
    //#pragma omp parallel for reduction(+:L)
    for (int i = 0; i < reflectance.rows; i++) {
        for (int j = 0; j < reflectance.cols; j++) {
            for (int k = 0; k < reflectance.depth; k++) {
                double val = abs(trueReflectance.getHostValue(i,j,k) - reflectance.getHostValue(i,j,k));
                L += val*val;
            }
        }
    }
    int N = reflectance.rows * reflectance.cols * reflectance.depth;
    return -0.5*beta*sqrt(L)/N;

}

void swapPointers(RealMatrix x, RealMatrix y) {
    double * temp = x.getHostPointer();
    x.setHostPointer(y.getHostPointer());
    y.setHostPointer(temp);
}

void computeReflectanceMatrix(Nanostructure proposedNanostructure, BioScat bioscat, RealMatrix reflectance, RealMatrix lambdas, RealMatrix betas, int num_lambdas, int num_betas, int n_obs)  {


    bioscat.getSegments(proposedNanostructure);

    
    for (int i = 0; i < num_lambdas; i++) {
        
        bioscat.reset();
        double lambda = lambdas.getHostValue(i);
        bioscat.prepareForward(lambda);

        for (int polarisation = 1; polarisation <= 2; polarisation++) {
            bioscat.forwardSolver(polarisation);
            bioscat.computeScatteredSubFields();
            bioscat.computeReflectedSubFields();
            bioscat.computeIncidentSubFields();
        }

        
        for (int j = 0; j < num_betas; j++) {
            double beta = betas.getHostValue(j);

            bioscat.computeScatteredFields(beta);

            bioscat.computeReflectedFields(beta);
        
            bioscat.computeIncidentFields(beta);
            
            bioscat.computeReflectance();

            for (int k = 0; k < n_obs; k++) {
                reflectance.setHostValue(k,i,j,bioscat.reflectance.getHostValue(k));
            }
        }
    }
}

void inverse(char * protein_structure, int num_segments, int total_grid_points, double * hyper, int num, int type_covfunc) {

    double start, stop; // Time measurement
    double Lprev, Lstar, alpha, u;
    double shift = 3e-8;
    double * temp;

    // Get number of data points
    int num_obs_points, num_lambdas, num_betas;
    getNumPoints(&num_obs_points, &num_lambdas, &num_betas, protein_structure, total_grid_points);
    
    // Allocate matrices for data
    RealMatrix trueReflectance = RealMatrix(num_obs_points, num_betas, num_lambdas);
    RealMatrix lambdas         = RealMatrix(num_lambdas);
    RealMatrix betas           = RealMatrix(num_betas);
    RealMatrix x_obs           = RealMatrix(num_obs_points);
    RealMatrix y_obs           = RealMatrix(num_obs_points);

    // Load data into matrices
    loadData(trueReflectance, lambdas, betas, x_obs, y_obs, protein_structure, total_grid_points);
   
    // Initialize some of the proposed nanostructures
    Nanostructure proposedNanostructure = Nanostructure(total_grid_points);
    initCoordInNanostructure(proposedNanostructure, protein_structure, total_grid_points);
    cudaDeviceSynchronize();
    // Setup Gaussian Process for realisations
    start = omp_get_wtime();
    GaussianProcess GP = GaussianProcess(total_grid_points, hyper, num, type_covfunc);
    stop = omp_get_wtime();

    printf("Initialization and allocation: %.4f seconds\n\n", stop - start);

    start = omp_get_wtime();
    GP.covariance_matrix();
    cudaDeviceSynchronize();
    stop = omp_get_wtime();

    printf("Computing covariance matrix: %.4f seconds\n\n", stop - start);

    start = omp_get_wtime();
    GP.cholesky();
    cudaDeviceSynchronize();
    stop = omp_get_wtime();

    printf("Cholesky factorization: %.4f seconds\n\n", stop - start);

    // Seed the random number generator with the current time
    srand(time(NULL));

    RealMatrix f            = RealMatrix(total_grid_points);
    RealMatrix fstar        = RealMatrix(total_grid_points);
    RealMatrix phi          = RealMatrix(total_grid_points);
    RealMatrix reflectance  = RealMatrix(num_obs_points, num_betas, num_lambdas);
    printf("HEJ %d %d %d\n",num_obs_points,num_betas,num_lambdas);
    
    //GP.realisation();

    for (int i = 0; i < total_grid_points; i++) {
        f.setHostValue(i,0.0);
    }

    for (int i = 0; i < total_grid_points; i++) {
        proposedNanostructure.f.setHostValue(i,f.getHostValue(i) + shift);
    }

    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points);

    //x_obs.print();
    //y_obs.print();


    bioscat.setupObservationPoints(x_obs.getHostPointer(), y_obs.getHostPointer(), num_obs_points);
    
    int acc =0;
    computeReflectanceMatrix(proposedNanostructure, bioscat, reflectance, lambdas, betas, num_lambdas, num_betas, num_obs_points);
    Lprev = computeLogLikelihood(trueReflectance, reflectance);
    printf("Lprev = %e\n", Lprev);
    computeReflectanceMatrix(proposedNanostructure, bioscat, reflectance, lambdas, betas, num_lambdas, num_betas, num_obs_points);
    Lprev = computeLogLikelihood(trueReflectance, reflectance);
    printf("Lprev = %e\n", Lprev);

    FILE *file;
    file = fopen("../../../Results/inverse/output.txt", "w");
    if (file == NULL) {
        perror("Error opening output file");
        return;
    }

    for (int j = 0; j < f.rows; j++) {
        fprintf(file, "%e ", f.getHostValue(j));
    }
    fprintf(file, "\n");

    double delta = 0.01;
    for (int n = 0; n < 200; n++) {
        GP.realisation();

        for (int i = 0; i < total_grid_points; i++) {
            phi.setHostValue(i,GP.p_h[i]*1e-8);
        }

        for (int i = 0; i < total_grid_points; i++) {
            fstar.setHostValue(i,sqrt(1 - 2*delta)*f.getHostValue(i) + sqrt(2*delta)*phi.getHostValue(i));
        }
        for (int i = 0; i < total_grid_points; i++) {
            proposedNanostructure.f.setHostValue(i,fstar.getHostValue(i) + shift);
        }

        computeReflectanceMatrix(proposedNanostructure, bioscat, reflectance, lambdas, betas, num_lambdas, num_betas, num_obs_points);

        // Compute log-likelihood
        Lstar = computeLogLikelihood(trueReflectance, reflectance);
        printf("(Lstar, Lprev) = (%f, %f)\n", Lstar, Lprev);

        // Compute probability of accepting curve
        alpha = std::min((double)1,exp(Lstar-Lprev));
        printf("alpha = %f\n", alpha);

        // Generate random number
        u = ((double) rand())/((double) RAND_MAX);
        printf("u = %f\n",u);

        if (u < alpha) {
            printf("ACCEPTED %d\n",acc++);
            //swapPointers(f, fstar);
            temp = f.getHostPointer();
            f.setHostPointer(fstar.getHostPointer());
            fstar.setHostPointer(temp);
            Lprev = Lstar;
            for (int j = 0; j < f.rows; j++) {
                fprintf(file, "%e ", f.getHostValue(j));
            }
            fprintf(file, "\n");
        }   
    }

    fclose(file);

    trueReflectance.free();
    lambdas.free();
    betas.free();
    x_obs.free();
    y_obs.free();
    GP.free();
    //bioscat.free();




}


}