#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <omp.h>
//#include <cblas.h>
#include <math.h>
#include "../../lib/BioScat.h"
#include "../../lib/Segment.h"
#include "../../lib/utils/RealMatrix.h"
extern "C" {
using namespace std;


__global__ void computeScatteredSubFieldsKernelCombine(ComplexMatrix field, ComplexMatrix subField, int rows) {
    int k = threadIdx.x + blockIdx.x * blockDim.x;
    if (k < rows) {
        field.setDeviceRealValue(k, field.getDeviceRealValue(k) + subField.getDeviceRealValue(k));
        field.setDeviceImagValue(k, field.getDeviceImagValue(k) + subField.getDeviceImagValue(k));
    }
}


void BioScat::computeScatteredSubFields() {

    int n = x_obs.rows;

    for (int i = 0; i < num_segments; i++) {
        
        double start_inner = omp_get_wtime();
        
        segments[i].computeScatteredFieldMatrices(x_obs, y_obs);

        segments[i].computeScatteredFields();

        double end_inner = omp_get_wtime();

        if (printOutput) printf("\nIt took %.4e seconds to compute the scattered field matrices for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }

    double start = omp_get_wtime();
    
    if (deviceComputation) {
     
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((n + dimBlock.x - 1)/dimBlock.x);

        if (polarisation == 1) {
        
            for (int i = 0; i < num_segments; i++) {
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(E_scat_pol[0].z, segments[i].E_scat.z, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(H_scat_pol[0].x, segments[i].H_scat.x, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(H_scat_pol[0].y, segments[i].H_scat.y, n);
                cudaDeviceSynchronize();
            }
           
        }
        else if (polarisation == 2) {

            for (int i = 0; i < num_segments; i++) {
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(H_scat_pol[1].z, segments[i].H_scat.z, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(E_scat_pol[1].x, segments[i].E_scat.x, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(E_scat_pol[1].y, segments[i].E_scat.y, n);
                cudaDeviceSynchronize();
            }

        }
        

    }
    else {
        if (polarisation == 1) {
            for (int k = 0; k < n; k++) {
                double Ez_real, Ez_imag, Hx_real, Hx_imag, Hy_real, Hy_imag, C_real, C_imag;
                Ez_real = 0.0; Hx_real = 0.0; Hy_real = 0.0;
                Ez_imag = 0.0; Hx_imag = 0.0; Hy_imag = 0.0;
                for (int i = 0; i < num_segments; i++) {
                    // Get real values
                    Ez_real += segments[i].E_scat.z.getHostRealValue(k);
                    Hx_real += segments[i].H_scat.x.getHostRealValue(k);
                    Hy_real += segments[i].H_scat.y.getHostRealValue(k);

                    // Get imagninary values
                    Ez_imag += segments[i].E_scat.z.getHostImagValue(k);
                    Hx_imag += segments[i].H_scat.x.getHostImagValue(k);
                    Hy_imag += segments[i].H_scat.y.getHostImagValue(k);
                    
                }

                E_scat_pol[0].z.setHostRealValue(k, Ez_real);
                H_scat_pol[0].x.setHostRealValue(k, Hx_real);
                H_scat_pol[0].y.setHostRealValue(k, Hy_real);
                E_scat_pol[0].z.setHostImagValue(k, Ez_imag);
                H_scat_pol[0].x.setHostImagValue(k, Hx_imag);
                H_scat_pol[0].y.setHostImagValue(k, Hy_imag);
            }
            
        }
        else if (polarisation == 2) {
            for (int k = 0; k < n; k++) {
                double Hz_real, Hz_imag, Ex_real, Ex_imag, Ey_real, Ey_imag, C_real, C_imag;
                Hz_real = 0.0; Ex_real = 0.0; Ey_real = 0.0;
                Hz_imag = 0.0; Ex_imag = 0.0; Ey_imag = 0.0;
                for (int i = 0; i < num_segments; i++) {

                        // Computing real values
                        Hz_real += segments[i].H_scat.z.getHostRealValue(k);
                        Ex_real += segments[i].E_scat.x.getHostRealValue(k);
                        Ey_real += segments[i].E_scat.y.getHostRealValue(k);

                        // Computing complex values
                        Hz_imag += segments[i].H_scat.z.getHostImagValue(k);
                        Ex_imag += segments[i].E_scat.x.getHostImagValue(k);
                        Ey_imag += segments[i].E_scat.y.getHostImagValue(k);
                
                }

                H_scat_pol[1].z.setHostRealValue(k, Hz_real);
                E_scat_pol[1].x.setHostRealValue(k, Ex_real);
                E_scat_pol[1].y.setHostRealValue(k, Ey_real);
                H_scat_pol[1].z.setHostImagValue(k, Hz_imag);
                E_scat_pol[1].x.setHostImagValue(k, Ex_imag);
                E_scat_pol[1].y.setHostImagValue(k, Ey_imag);
            }  
        }
    }

    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %.4e seconds to compute the combined scattered fields in the observation points.\n\n",end-start);
    for (int i = 0; i < num_segments; i++) {
        //segments[i].freeScatteredSubFields();
    }

}

void BioScat::computeInteriorSubFields() {

    int n = x_obs.rows;

    for (int i = 0; i < num_segments; i++) {
        
        double start_inner = omp_get_wtime();
        
        segments[i].computeInteriorFieldMatrices(x_obs, y_obs);

        segments[i].computeInteriorFields();

        double end_inner = omp_get_wtime();

        if (printOutput) printf("\nIt took %.4e seconds to compute the interior field matrices for segment %d.\n\n",end_inner - start_inner, i + 1);
        
    }

    double start = omp_get_wtime();
    
    if (deviceComputation) {
     
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((n + dimBlock.x - 1)/dimBlock.x);

        if (polarisation == 1) {
        
            for (int i = 0; i < num_segments; i++) {
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(E_int_pol[0].z, segments[i].E_int.z, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(H_int_pol[0].x, segments[i].H_int.x, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(H_int_pol[0].y, segments[i].H_int.y, n);
                cudaDeviceSynchronize();
            }
           
        }
        else if (polarisation == 2) {

            for (int i = 0; i < num_segments; i++) {
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(H_int_pol[1].z, segments[i].H_int.z, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(E_int_pol[1].x, segments[i].E_int.x, n);
                computeScatteredSubFieldsKernelCombine<<<dimGrid, dimBlock>>>(E_int_pol[1].y, segments[i].E_int.y, n);
                cudaDeviceSynchronize();
            }

        }
        

    }
    else {
        if (polarisation == 1) {
            for (int k = 0; k < n; k++) {
                double Ez_real, Ez_imag, Hx_real, Hx_imag, Hy_real, Hy_imag, C_real, C_imag;
                Ez_real = 0.0; Hx_real = 0.0; Hy_real = 0.0;
                Ez_imag = 0.0; Hx_imag = 0.0; Hy_imag = 0.0;
                for (int i = 0; i < num_segments; i++) {
                    // Get real values
                    Ez_real += segments[i].E_int.z.getHostRealValue(k);
                    Hx_real += segments[i].H_int.x.getHostRealValue(k);
                    Hy_real += segments[i].H_int.y.getHostRealValue(k);

                    // Get imagninary values
                    Ez_imag += segments[i].E_int.z.getHostImagValue(k);
                    Hx_imag += segments[i].H_int.x.getHostImagValue(k);
                    Hy_imag += segments[i].H_int.y.getHostImagValue(k);
                    
                }

                E_int_pol[0].z.setHostRealValue(k, Ez_real);
                H_int_pol[0].x.setHostRealValue(k, Hx_real);
                H_int_pol[0].y.setHostRealValue(k, Hy_real);
                E_int_pol[0].z.setHostImagValue(k, Ez_imag);
                H_int_pol[0].x.setHostImagValue(k, Hx_imag);
                H_int_pol[0].y.setHostImagValue(k, Hy_imag);
            }
            
        }
        else if (polarisation == 2) {
            for (int k = 0; k < n; k++) {
                double Hz_real, Hz_imag, Ex_real, Ex_imag, Ey_real, Ey_imag, C_real, C_imag;
                Hz_real = 0.0; Ex_real = 0.0; Ey_real = 0.0;
                Hz_imag = 0.0; Ex_imag = 0.0; Ey_imag = 0.0;
                for (int i = 0; i < num_segments; i++) {

                        // Computing real values
                        Hz_real += segments[i].H_int.z.getHostRealValue(k);
                        Ex_real += segments[i].E_int.x.getHostRealValue(k);
                        Ey_real += segments[i].E_int.y.getHostRealValue(k);

                        // Computing complex values
                        Hz_imag += segments[i].H_int.z.getHostImagValue(k);
                        Ex_imag += segments[i].E_int.x.getHostImagValue(k);
                        Ey_imag += segments[i].E_int.y.getHostImagValue(k);
                
                }

                H_int_pol[1].z.setHostRealValue(k, Hz_real);
                E_int_pol[1].x.setHostRealValue(k, Ex_real);
                E_int_pol[1].y.setHostRealValue(k, Ey_real);
                H_int_pol[1].z.setHostImagValue(k, Hz_imag);
                E_int_pol[1].x.setHostImagValue(k, Ex_imag);
                E_int_pol[1].y.setHostImagValue(k, Ey_imag);
            }  
        }
    }

    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %.4e seconds to compute the combined interior fields in the observation points.\n\n",end-start);
    for (int i = 0; i < num_segments; i++) {
        //segments[i].freeScatteredSubFields();
    }

}


__global__ void copyRealKernel(ComplexMatrix output, ComplexMatrix input, int n) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n) {
        output.setDeviceRealValue(j,input.getDeviceRealValue(j));
    }
}

__global__ void copyImagKernel(ComplexMatrix output, ComplexMatrix input, int n) {
    int j = threadIdx.x + blockIdx.x * blockDim.x;
    if (j < n) {
        output.setDeviceImagValue(j,input.getDeviceImagValue(j));
    }
}


void BioScat::computeIncidentSubFields() {
    int n = x_obs.rows;
    double val;
    double start = omp_get_wtime();
    segments[0].computeIncidentFields(y_obs);
    if (deviceComputation) {
        // Blocks and threads
        dim3 dimBlock(256);
        dim3 dimGrid((n + dimBlock.x - 1)/dimBlock.x);
        if (polarisation == 1) {
            copyRealKernel<<<dimGrid, dimBlock>>>(E_inc_pol[0].z, segments[0].E_inc.z, n);
            copyImagKernel<<<dimGrid, dimBlock>>>(E_inc_pol[0].z, segments[0].E_inc.z, n);
            copyRealKernel<<<dimGrid, dimBlock>>>(H_inc_pol[0].x, segments[0].H_inc.x, n);
            copyImagKernel<<<dimGrid, dimBlock>>>(H_inc_pol[0].x, segments[0].H_inc.x, n);
        }
        else if (polarisation == 2) {
            copyRealKernel<<<dimGrid, dimBlock>>>(H_inc_pol[1].z, segments[0].H_inc.z, n);
            copyImagKernel<<<dimGrid, dimBlock>>>(H_inc_pol[1].z, segments[0].H_inc.z, n);
            copyRealKernel<<<dimGrid, dimBlock>>>(E_inc_pol[1].x, segments[0].E_inc.x, n);
            copyImagKernel<<<dimGrid, dimBlock>>>(E_inc_pol[1].x, segments[0].E_inc.x, n);
        }

        cudaDeviceSynchronize();
    }
    else {
        if (polarisation == 1) {
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) E_inc_pol[0].z.setHostRealValue(k, segments[0].E_inc.z.getHostRealValue(k));
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) E_inc_pol[0].z.setHostImagValue(k, segments[0].E_inc.z.getHostImagValue(k));
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) H_inc_pol[0].x.setHostRealValue(k, segments[0].H_inc.x.getHostRealValue(k));
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) H_inc_pol[0].x.setHostImagValue(k, segments[0].H_inc.x.getHostImagValue(k));
        }
        else if (polarisation == 2) {
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) H_inc_pol[1].z.setHostRealValue(k, segments[0].H_inc.z.getHostRealValue(k));
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) H_inc_pol[1].z.setHostImagValue(k, segments[0].H_inc.z.getHostImagValue(k));
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) E_inc_pol[1].x.setHostRealValue(k, segments[0].E_inc.x.getHostRealValue(k));
            //#pragma omp parallel for
            for (int k = 0; k < n; k++) E_inc_pol[1].x.setHostImagValue(k, segments[0].E_inc.x.getHostImagValue(k));
        }
    }
    double end = omp_get_wtime();
    if (printOutput) printf("\nIt took %.4e seconds to compute the combined incident fields in the observation points.\n\n",end-start);
    
}


}