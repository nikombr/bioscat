
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
using namespace std;

BioScat::BioScat(char* protein_structure, int num_segments, int total_grid_points) {

    this->protein_structure = protein_structure;
    this->num_segments = num_segments;
    this->total_grid_points = total_grid_points;

}

BioScat::~BioScat() {

    for (int i = 0; i < num_segments; i++) {
        segments[i].free();
    }

    delete[] segments;
    reflectance.free();

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
    for (int i = 0; i < num_segments; i++) {
        segments[i].polarisation = polarisation;

        start_inner = omp_get_wtime();
        segments[i].computeFieldsForLinearSystem();
        segments[i].setupRightHandSide();
        printf("HEJ2\n");
        segments[i].setupSystemMatrix();
        segments[i].freeScatteredFields();
        segments[i].freeInteriorFields();
        segments[i].freeIncidentFields();
        segments[i].freeReflectedFields();
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
        // Initialize all arrays to zero on the host
        E_scat_pol[i].setHostZero();
        H_scat_pol[i].setHostZero();
        E_inc_pol[i].setHostZero();
        H_inc_pol[i].setHostZero();
        E_ref_pol[i].setHostZero();
        H_ref_pol[i].setHostZero();
    }

    // Allocate reflectance array
    reflectance = RealMatrix(n);
}

void BioScat::computeScatteredFields() {

    computeScatteredFields(beta);
}

void BioScat::computeScatteredSubFields() {

    int n = x_obs.rows;

    double start = omp_get_wtime();

    for (int i = 0; i < num_segments; i++) {
        
        double start_inner = omp_get_wtime();
        
        segments[i].computeScatteredFieldMatrices(x_obs, y_obs, false);

        double end_inner = omp_get_wtime();
        printf("\nIt took %.4e seconds to compute the scattered field matrices for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }
    
    if (polarisation == 1) {
        //#pragma omp parallel for
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

            E_scat_pol[0].z.setHostRealValue(k, val1);
            H_scat_pol[0].x.setHostRealValue(k, val2);
            H_scat_pol[0].y.setHostRealValue(k, val3);
        }
    
        //#pragma omp parallel for
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
            E_scat_pol[0].z.setHostImagValue(k, val1);
            H_scat_pol[0].x.setHostImagValue(k, val2);
            H_scat_pol[0].y.setHostImagValue(k, val3);
            
        }
           
    }
    else if (polarisation == 2) {
        //#pragma omp parallel for
        for (int k = 0; k < n; k++) {
            double val1, val2, val3, C_real, C_imag;
            val1 = 0.0;
            val2 = 0.0;
            val3 = 0.0;
            for (int i = 0; i < num_segments; i++) {
                
                for (int j = 0; j < segments[i].n_int; j++) {
                    C_real = segments[i].C.getHostRealValue(j);
                    C_imag = segments[i].C.getHostImagValue(j);
                    val1 += segments[i].H_scat_matrix.z.getHostRealValue(k,j) * C_real;
                    val1 -= segments[i].H_scat_matrix.z.getHostImagValue(k,j) * C_imag;
                    val2 += segments[i].E_scat_matrix.x.getHostRealValue(k,j) * C_real;
                    val2 -= segments[i].E_scat_matrix.x.getHostImagValue(k,j) * C_imag;
                    val3 += segments[i].E_scat_matrix.y.getHostRealValue(k,j) * C_real;
                    val3 -= segments[i].E_scat_matrix.y.getHostImagValue(k,j) * C_imag;
            
                }
                
            }

            H_scat_pol[1].z.setHostRealValue(k, val1);
            E_scat_pol[1].x.setHostRealValue(k, val2);
            E_scat_pol[1].y.setHostRealValue(k, val3);
        }
    
        //#pragma omp parallel for
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
                    val1 += segments[i].H_scat_matrix.z.getHostRealValue(k,j) * C_imag;
                    val1 += segments[i].H_scat_matrix.z.getHostImagValue(k,j) * C_real;
                    val2 += segments[i].E_scat_matrix.x.getHostRealValue(k,j) * C_imag;
                    val2 += segments[i].E_scat_matrix.x.getHostImagValue(k,j) * C_real;
                    val3 += segments[i].E_scat_matrix.y.getHostRealValue(k,j) * C_imag;
                    val3 += segments[i].E_scat_matrix.y.getHostImagValue(k,j) * C_real;
                }
                //}
            }
            H_scat_pol[1].z.setHostImagValue(k, val1);
            E_scat_pol[1].x.setHostImagValue(k, val2);
            E_scat_pol[1].y.setHostImagValue(k, val3);
            
        }
           
    }

    double end = omp_get_wtime();
    printf("\nIt took %.4e seconds to compute the combined scattered fields in the observation points.\n\n",end-start);
    for (int i = 0; i < num_segments; i++) {
        segments[i].freeScatteredFields();
    }

}

void combinePolarisation(Field * pol, Field combined, double beta) {
    double cosBeta = cos(beta);
    double sinBeta = sin(beta);
    for (int i = 0; i < combined.x.rows; i++) {
        combined.x.setHostRealValue(i, pol[0].x.getHostRealValue(i)*cosBeta + pol[1].x.getHostRealValue(i)*sinBeta);
        combined.y.setHostRealValue(i, pol[0].y.getHostRealValue(i)*cosBeta + pol[1].y.getHostRealValue(i)*sinBeta);
        combined.z.setHostRealValue(i, pol[0].z.getHostRealValue(i)*cosBeta + pol[1].z.getHostRealValue(i)*sinBeta);
        combined.x.setHostImagValue(i, pol[0].x.getHostImagValue(i)*cosBeta + pol[1].x.getHostImagValue(i)*sinBeta);
        combined.y.setHostImagValue(i, pol[0].y.getHostImagValue(i)*cosBeta + pol[1].y.getHostImagValue(i)*sinBeta);
        combined.z.setHostImagValue(i, pol[0].z.getHostImagValue(i)*cosBeta + pol[1].z.getHostImagValue(i)*sinBeta);
    }

}

void BioScat::computeScatteredFields(double beta) {
    combinePolarisation(E_scat_pol, E_scat, beta);
    combinePolarisation(H_scat_pol, H_scat, beta);
    
}

void BioScat::computeIncidentFields() {
    
    computeIncidentFields(beta);
}

void BioScat::computeIncidentSubFields() {
    int n = x_obs.rows;
    double val;
    double start = omp_get_wtime();
    segments[0].computeIncidentFieldVectors(y_obs);
    
    if (polarisation == 1) {
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.z.getHostRealValue(k);
            E_inc_pol[0].z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.z.getHostImagValue(k);
            E_inc_pol[0].z.setHostImagValue(k, val);
        }

        for (int k = 0; k < n; k++) {
            val = segments[0].H_inc_vector.x.getHostRealValue(k);
            H_inc_pol[0].x.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].H_inc_vector.x.getHostImagValue(k);
            H_inc_pol[0].x.setHostImagValue(k, val);
        }
    }
    else if (polarisation == 2) {
        for (int k = 0; k < n; k++) {
            val = segments[0].H_inc_vector.z.getHostRealValue(k);
            H_inc_pol[1].z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].H_inc_vector.z.getHostImagValue(k);
            H_inc_pol[1].z.setHostImagValue(k, val);
        }

        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.x.getHostRealValue(k);
            E_inc_pol[1].x.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_inc_vector.x.getHostImagValue(k);
            E_inc_pol[1].x.setHostImagValue(k, val);
        }
    }
    double end = omp_get_wtime();
    printf("\nIt took %.4e seconds to compute the combined incident fields in the observation points.\n\n",end-start);
    
    segments[0].freeIncidentFields();
    
}

void BioScat::computeIncidentFields(double beta) {

    combinePolarisation(E_inc_pol, E_inc, beta);
    combinePolarisation(H_inc_pol, H_inc, beta);

    
}

void BioScat::computeReflectedFields() {
    computeReflectedFields(beta);
}

void BioScat::computeReflectedSubFields() {

    int n = x_obs.rows;
    double val;
    double start = omp_get_wtime();
    segments[0].computeReflectedFieldVectors(y_obs);

    if (polarisation == 1) {
        for (int k = 0; k < n; k++) {
            val = segments[0].E_ref_vector.z.getHostRealValue(k);
            E_ref_pol[0].z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_ref_vector.z.getHostImagValue(k);
            E_ref_pol[0].z.setHostImagValue(k, val);
        }

        for (int k = 0; k < n; k++) {
            val = segments[0].H_ref_vector.x.getHostRealValue(k);
            H_ref_pol[0].x.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].H_ref_vector.x.getHostImagValue(k);
            H_ref_pol[0].x.setHostImagValue(k, val);
        }
    }
    else if (polarisation == 2) {
        for (int k = 0; k < n; k++) {
            val = segments[0].H_ref_vector.z.getHostRealValue(k);
            H_ref_pol[1].z.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].H_ref_vector.z.getHostImagValue(k);
            H_ref_pol[1].z.setHostImagValue(k, val);
        }

        for (int k = 0; k < n; k++) {
            val = segments[0].E_ref_vector.x.getHostRealValue(k);
            E_ref_pol[1].x.setHostRealValue(k, val);
        }
        for (int k = 0; k < n; k++) {
            val = segments[0].E_ref_vector.x.getHostImagValue(k);
            E_ref_pol[1].x.setHostImagValue(k, val);
        }
    }
    double end = omp_get_wtime();
    printf("\nIt took %.4e seconds to compute the combined reflected fields in the observation points.\n\n",end-start);
    
    segments[0].freeReflectedFields();
    
}

void BioScat::computeReflectedFields(double beta) {

    combinePolarisation(E_ref_pol, E_ref, beta);
    combinePolarisation(H_ref_pol, H_ref, beta);

    
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

}