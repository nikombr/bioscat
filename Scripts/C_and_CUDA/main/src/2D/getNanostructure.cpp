#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include <fstream>
extern "C" {
#include "../../lib/BioScat.h"
using namespace std;


void BioScat::getNanostructure() {

    char filename[256];
    Nanostructure nanostructure(total_grid_points); 

    sprintf(filename, "../../../../Data/nanostructures/%s_2D_x_%d.txt", protein_structure, total_grid_points);

    FILE *file;
    file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }
    for (int i = 0; i < total_grid_points; i++) {
        fprintf(file, "%lf,", nanostructure.x.val_h[i]);  // Reading each value into the array
    }

    sprintf(filename, "../../../../Data/nanostructures/%s_2D_f_%d.txt", protein_structure, total_grid_points);

    file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        return;
    }

    for (int i = 0; i < total_grid_points; i++) {
        fprintf(file, "%lf,", nanostructure.f.val_h[i]);  // Reading each value into the array
    }



}



}