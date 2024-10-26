
#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../lib/BioScat.h"
#include "../lib/2D/Segment.h"
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

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].setup(this->nanostructure, i, total_grid_points, num_segments);
    }

}

void BioScat::getSegments(Nanostructure nanostructure) {

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].setup(nanostructure, i, total_grid_points, num_segments);
    }

}

void BioScat::forwardSolver(int scenario) {

    for (int i = 0; i < num_segments; i++) {
        segments[i].scenario = scenario;
        segments[i].computeFieldsForLinearSystem();
        segments[i].setupRightHandSide();
        segments[i].setupSystemSubMatrices();
        segments[i].solveLinearSystem();
        

        
    }

}

}