#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
#include <math.h>
#include <stdexcept>
#include "../../lib/utils/ComplexMatrix.h"
#include "../../lib/BioScat.h"
#include "../../lib/reflectance/computeReflectanceMatrix.h"
extern "C" {

using namespace std;
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

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

void getNumPoints(int * n_obs, int * num_lambdas, int * num_betas, char * protein_structure,char * datatype) {
    // Get number of data points
    char filename[256];
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");

    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/phi.txt",dir,datatype,protein_structure);

    * n_obs = countNumbersInFile(filename);

    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/lambdas.txt",dir,datatype,protein_structure);

    * num_lambdas = countNumbersInFile(filename);
   
    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/betas.txt",dir,datatype,protein_structure);

    * num_betas = countNumbersInFile(filename);
    
}

void loadData(RealMatrix trueReflectance, RealMatrix lambdas, RealMatrix betas, RealMatrix phi, char * protein_structure, char * datatype) {
    char filename[256];
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");
    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/reflectance.txt",dir,datatype,protein_structure);
    trueReflectance.loadVector(filename);
    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/lambdas.txt",dir,datatype,protein_structure);
    lambdas.loadVector(filename);
    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/betas.txt",dir,datatype,protein_structure);
    betas.loadVector(filename);
    sprintf(filename, "%s/Data/artificial_data/comsol/%s/%s/phi.txt",dir,datatype,protein_structure);
    phi.loadVector(filename);

}

void initCoordInNanostructure(Nanostructure nanostructure, char * protein_structure, int total_grid_points) {
    double val;
    char filename[256];
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");
    sprintf(filename, "%s/Data/nanostructures/2D/%s_x_%d.txt", dir, protein_structure, total_grid_points);

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

    nanostructure.x.toDevice();
    cudaDeviceSynchronize();
}

double computeLogLikelihood(RealMatrix trueReflectance, RealMatrix reflectance, bool deviceComputation, double gamma) {

    if (deviceComputation) {
        reflectance.toHost();
    }

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
    //int N = reflectance.rows * reflectance.cols * reflectance.depth;
    return -0.5*gamma*sqrt(L);
    
    

}

void swapPointers(RealMatrix x, RealMatrix y) {
    double * temp = x.getHostPointer();
    x.setHostPointer(y.getHostPointer());
    y.setHostPointer(temp);
}


double computeInverseStep(Nanostructure proposedNanostructure, BioScat bioscat, RealMatrix reflectance, RealMatrix trueReflectance, RealMatrix lambdas, RealMatrix betas, int num_lambdas, int num_betas, int n_obs, bool deviceComputation, double gamma) {
    computeReflectanceMatrix(proposedNanostructure, bioscat, reflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation);
    //printf("test01 %f\n",computeLogLikelihood(trueReflectance, reflectance, deviceComputation));
    //computeReflectanceMatrix(proposedNanostructure, bioscat, reflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation);
    //printf("test02 %f\n",computeLogLikelihood(trueReflectance, reflectance, deviceComputation));
    //computeReflectanceMatrix(proposedNanostructure, bioscat, reflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation);
    //printf("test03 %f\n",computeLogLikelihood(trueReflectance, reflectance, deviceComputation));
    return computeLogLikelihood(trueReflectance, reflectance, deviceComputation, gamma);
}

void inverse(char * protein_structure, int num_segments, int total_grid_points, double * hyper, int num, int type_covfunc, double delta, int maxiter, char * filename, double decay_rate, double gamma, int fine_tuning_int, char * datatype, int chainNum) {

    double start, stop; // Time measurement
    double Lprev, Lstar, alpha, u, logPrior;
    double shift = 1.5e-6;
    double scale = 1e-8;
    double * temp_h, * temp_d;
    bool deviceComputation = true;
    double deltastart = delta;
    bool fine_tuning = fine_tuning_int == 1 ? true : false;
    
    // Files for output
    FILE *file, *logfile, *logfile_accepted;
    char current_file_name[256];
    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");
    printf("%s\n",datatype);
    printf("%s\n",filename);
    sprintf(current_file_name,"%s/Results/inverse/%s/%s_output.txt",dir,datatype,filename);
    file = fopen(current_file_name, "w");
    if (file == NULL) {
        printf("Error opening output file 1: %s\n",current_file_name);
        return;
    }
    sprintf(current_file_name,"%s/Results/inverse/%s/%s_log.txt",dir,datatype,filename);
    logfile = fopen(current_file_name, "w");
    if (file == NULL) {
        printf("Error opening output file 2: %s\n",current_file_name);
        return;
    }
    sprintf(current_file_name,"%s/Results/inverse/%s/%s_log_accepted.txt",dir,datatype,filename);
    logfile_accepted = fopen(current_file_name, "w");
    if (file == NULL) {
        printf("Error opening output file 3: %s\n",current_file_name);
        return;
    }

    // Get number of data points
    int n_obs, num_lambdas, num_betas;
    getNumPoints(&n_obs, &num_lambdas, &num_betas, protein_structure, datatype);
    printf("(n_obs, n_lambdas, n_beta) = (%d, %d, %d)\n",n_obs,num_lambdas,num_betas);

    // Allocate matrices for data
    RealMatrix trueReflectance = RealMatrix(n_obs, num_betas, num_lambdas);
    RealMatrix lambdas         = RealMatrix(num_lambdas);
    RealMatrix betas           = RealMatrix(num_betas);
    RealMatrix angles          = RealMatrix(n_obs);
    RealMatrix f               = RealMatrix(total_grid_points);
    RealMatrix fstar           = RealMatrix(total_grid_points);
    RealMatrix phi             = RealMatrix(total_grid_points);
    RealMatrix reflectance     = RealMatrix(n_obs, num_betas, num_lambdas);

    // Load data into matrices
    loadData(trueReflectance, lambdas, betas, angles, protein_structure, datatype);
   
    // Initialize some of the proposed nanostructures
    Nanostructure proposedNanostructure = Nanostructure(total_grid_points);
    initCoordInNanostructure(proposedNanostructure, protein_structure, total_grid_points);
    
    // Setup Gaussian Process for realisations
    start = omp_get_wtime();
    GaussianProcess GP = GaussianProcess(total_grid_points, hyper, num, type_covfunc);
    stop = omp_get_wtime();

    printf("Initialization and allocation: %.4f seconds\n\n", stop - start);

    start = omp_get_wtime();
    GP.covariance_matrix();
    stop = omp_get_wtime();

    printf("Computing covariance matrix: %.4f seconds\n\n", stop - start);

    start = omp_get_wtime();
    GP.compute_inverse();
    stop = omp_get_wtime();

    printf("Computing inverse of covariance matrix: %.4f seconds\n\n", stop - start);

    start = omp_get_wtime();
    GP.cholesky();
    stop = omp_get_wtime();

    printf("Cholesky factorization: %.4f seconds\n\n", stop - start);

    
    // Seed the random number generator with the current time
    srand(time(NULL));
    
    GP.realisation();
    // Send to host
    cudaMemcpy(GP.p_h, GP.p_d, total_grid_points * sizeof(double), cudaMemcpyDeviceToHost); // remove at some point perhaps
    for (int i = 0; i < total_grid_points; i++) {
        f.setHostValue(i,GP.p_h[i]);
    }

    /*for (int i = 0; i < total_grid_points; i++) {
        f.setHostValue(i,0.0);
    }*/

    for (int i = 0; i < total_grid_points; i++) {
        proposedNanostructure.f.setHostValue(i,f.getHostValue(i)*scale + shift);
    }
    bool printOutput = false;
    BioScat bioscat = BioScat(protein_structure, num_segments, total_grid_points, deviceComputation, angles.getHostPointer(), n_obs, printOutput);

    int acc = 0; // Count how many curves are accepted

    f.toDevice();
    logPrior = GP.compute_prior(f.getDevicePointer());

    Lprev = computeInverseStep(proposedNanostructure, bioscat, reflectance, trueReflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation, gamma);
 
    for (int j = 0; j < f.rows; j++) {
        fprintf(file, "%e ", proposedNanostructure.f.getHostValue(j));
    }
    fprintf(file, "\n");
    fflush(file);
    
    fprintf(logfile, "accepted delta gamma Lprev Lstar alpha time\n");
    fprintf(logfile_accepted, "%e %e %e %e %e %e\n",Lprev, exp(Lprev),logPrior, exp(logPrior), Lprev + logPrior, exp(Lprev + logPrior));

    int status;
    for (int n = 0; n < maxiter; n++) {

        // Decay rate applied
        delta = std::max(deltastart/(1 + decay_rate*n), 1e-6);

        // Fine tuning
        if (fine_tuning && -Lprev < 25) {
            gamma = std::min(1.1*gamma,1e5);
            Lprev = computeInverseStep(proposedNanostructure, bioscat, reflectance, trueReflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation, gamma);
        }
    
        GP.realisation();
        // Send to host
        cudaMemcpy(GP.p_h, GP.p_d, total_grid_points * sizeof(double), cudaMemcpyDeviceToHost); // remove at some point perhaps
    

        for (int i = 0; i < total_grid_points; i++) {
            phi.setHostValue(i,GP.p_h[i]);
        }

        for (int i = 0; i < total_grid_points; i++) {
            fstar.setHostValue(i,sqrt(1 - 2*delta)*f.getHostValue(i) + sqrt(2*delta)*phi.getHostValue(i));
        }
        /*for (int i = 0; i < total_grid_points; i++) {
            proposedNanostructure.f.setHostValue(i,f.getHostValue(i)*scale + shift);
        }
        double Lhmm = computeInverseStep(proposedNanostructure, bioscat, reflectance, trueReflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation);*/
        for (int i = 0; i < total_grid_points; i++) {
            proposedNanostructure.f.setHostValue(i,fstar.getHostValue(i)*scale + shift);
        }
      
        double start = omp_get_wtime();
        Lstar = computeInverseStep(proposedNanostructure, bioscat, reflectance, trueReflectance, lambdas, betas, num_lambdas, num_betas, n_obs, deviceComputation, gamma);
        double end = omp_get_wtime();
        //printf("L = (%f, %f)\n", Lhmm, Lstar);
        /*if (status == EXIT_SUCCESS) {
        // Compute log-likelihood
        Lstar = computeLogLikelihood(trueReflectance, reflectance);
        if (minimum < 1e-8) Lstar *= 100;
        //printf("(Lstar, Lprev) = (%f, %f)\n", Lstar, Lprev);*/

        // Compute probability of accepting curve   
        alpha = std::min(1.0, exp(Lstar-Lprev));

        // Generate random number
        u = ((double) rand())/((double) RAND_MAX);

        if (u < alpha) {
            acc++;
            //printf("ACCEPTED %d\n",acc++);
            //swapPointers(f, fstar);
            temp_h = f.getHostPointer();
            temp_d = f.getDevicePointer();
            f.setHostPointer(fstar.getHostPointer());
            f.setDevicePointer(fstar.getDevicePointer());
            fstar.setHostPointer(temp_h);
            fstar.setDevicePointer(temp_d);
            Lprev = Lstar;
            for (int j = 0; j < f.rows; j++) {
                fprintf(file, "%e ", proposedNanostructure.f.getHostValue(j));
            }
            fprintf(file, "\n");
            fflush(file);
            f.toDevice();
            logPrior = GP.compute_prior(f.getDevicePointer());
            fprintf(logfile_accepted, "%e %e %e %e %e %e\n",Lprev, exp(Lprev),logPrior, exp(logPrior), Lprev + logPrior, exp(Lprev + logPrior));
            fflush(logfile_accepted);
        }

        fprintf(logfile, "%d %e %e %f %f %f %e\n", acc, delta, gamma, Lprev, Lstar, alpha, end - start);
        fflush(logfile);

        if (acc > 5000) break;
          
    }

    fclose(file);
    fclose(logfile);

    trueReflectance.free();
    lambdas.free();
    betas.free();
    phi.free();
    GP.free();
    bioscat.free();




}


}