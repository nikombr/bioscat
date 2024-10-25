#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <fstream>
extern "C" {
#include "../../lib/BioScat.h"
#include "../../lib/RealMatrix.h"
using namespace std;


void BioScat::getNanostructure() {

    char filename[256];
    nanostructure = Nanostructure(total_grid_points); 

    sprintf(filename, "../../../Data/nanostructures/%s_2D_x_%d.txt", protein_structure, total_grid_points);
    //printf("filename = %s\n",filename);

    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file 1");
        return;
    }
    for (int i = 0; i < total_grid_points; i++) {
        fscanf(file, "%lf,", &nanostructure.x.val_h[i]);  // Reading each value into the array
    }
    fclose(file);

    sprintf(filename, "../../../Data/nanostructures/%s_2D_f_%d.txt", protein_structure, total_grid_points);

    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file 2");
        return;
    }

    for (int i = 0; i < total_grid_points; i++) {
        fscanf(file, "%lf,", &nanostructure.f.val_h[i]);  // Reading each value into the array
    }

    fclose(file);

    // Move data to the device
    nanostructure.x.toDevice();
    nanostructure.f.toDevice();


}



}