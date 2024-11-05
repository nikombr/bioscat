#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

void Segment::computeIncidentFieldVectors(RealMatrix y) {
    
    int rows = y.rows;

    bool Ex_bool = polarisation == 1 ? false : true;
    bool Ez_bool = polarisation == 1 ? true  : false;
    bool Hx_bool = polarisation == 1 ? true  : false;
    bool Hz_bool = polarisation == 1 ? false : true;
    bool Ey_bool = false;
    bool Hy_bool = false;

    E_inc_vector = Field(rows);
    H_inc_vector = Field(rows);

    if (polarisation == 1) {
        for (int j = 0; j < rows; j++) E_inc_vector.z.setHostRealValue(j,                     cos(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_inc_vector.z.setHostImagValue(j,                     sin(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.x.setHostRealValue(j, -1/constants.eta0 * cos(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.x.setHostImagValue(j, -1/constants.eta0 * sin(constants.k0 * y.getHostValue(j)));
    }
    else {
        for (int j = 0; j < rows; j++) E_inc_vector.x.setHostRealValue(j,                    cos(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_inc_vector.x.setHostImagValue(j,                    sin(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.z.setHostRealValue(j, 1/constants.eta0 * cos(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_inc_vector.z.setHostImagValue(j, 1/constants.eta0 * sin(constants.k0 * y.getHostValue(j)));
    }


}

void Segment::computeReflectedFieldVectors(RealMatrix y) {
    
    int rows = y.rows;

    bool Ex_bool = polarisation == 1 ? false : true;
    bool Ez_bool = polarisation == 1 ? true  : false;
    bool Hx_bool = polarisation == 1 ? true  : false;
    bool Hz_bool = polarisation == 1 ? false : true;
    bool Ey_bool = false;
    bool Hy_bool = false;

    E_ref_vector = Field(rows);
    H_ref_vector = Field(rows);

    if (polarisation == 1) {
        for (int j = 0; j < rows; j++) {
            E_ref_vector.z.setHostRealValue(j,                     constants.Gamma_ref * cos(constants.k0 * y.getHostValue(j)));
            E_ref_vector.z.setHostImagValue(j,                    -constants.Gamma_ref * sin(constants.k0 * y.getHostValue(j)));
            H_ref_vector.x.setHostRealValue(j,  1/constants.eta0 * constants.Gamma_ref * cos(constants.k0 * y.getHostValue(j)));
            H_ref_vector.x.setHostImagValue(j, -1/constants.eta0 * constants.Gamma_ref * sin(constants.k0 * y.getHostValue(j)));
        }
    }
    else {
        for (int j = 0; j < rows; j++) E_ref_vector.x.setHostRealValue(j,                     constants.Gamma_ref * cos(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) E_ref_vector.x.setHostImagValue(j,                    -constants.Gamma_ref * sin(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_ref_vector.z.setHostRealValue(j, -1/constants.eta0 * constants.Gamma_ref * cos(constants.k0 * y.getHostValue(j)));
        for (int j = 0; j < rows; j++) H_ref_vector.z.setHostImagValue(j,  1/constants.eta0 * constants.Gamma_ref * sin(constants.k0 * y.getHostValue(j)));
    }


}
  

}