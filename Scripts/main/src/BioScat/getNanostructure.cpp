#include <stdlib.h>
#include <stdio.h>
//#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <fstream>
#include "../../lib/BioScat.h"
//#include "../lib/RealMatrix.h"
extern "C" {
using namespace std;


void BioScat::getNanostructure() {
    // Loads information about nanostructure and transfers to device

    // Initialize variables
    char filename[256];
    double val;

    // Initialize size of matrices to store informatino about nanostructure
    nanostructure = Nanostructure(total_grid_points); 

    char dir[256];
    sprintf(dir,"../../../../../../../work3/s194146/bioscatdata");

    sprintf(filename, "%s/Data/nanostructures/2D/%s_x_%d.txt", dir, protein_structure, total_grid_points);

    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file %s\n",filename);
        return;
    }
    for (int i = 0; i < total_grid_points; i++) {
        fscanf(file, "%lf,", &val);  // Reading each value into the array
        nanostructure.x.setHostValue(i, val);
    }
    fclose(file);

    sprintf(filename, "%s/Data/nanostructures/2D/%s_f_%d.txt", dir, protein_structure, total_grid_points);

    file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file %s\n",filename);
        return;
    }

    for (int i = 0; i < total_grid_points; i++) {
        fscanf(file, "%lf,", &val);  // Reading each value into the array
        nanostructure.f.setHostValue(i, val);
    }

    fclose(file);

    if (deviceComputation) {
        // Move data to the device
        nanostructure.x.toDevice();
        nanostructure.f.toDevice();
    }

}



}