#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/2D/Segment.h"
#include "../../lib/BioScat.h"
#include "../../lib/RealMatrix.h"
using namespace std;


void BioScat::getSegments() {

    this->segments = new Segment[num_segments];

    for (int i = 0; i < num_segments; i++) {
        segments[i].setup(nanostructure, i, total_grid_points, num_segments);
    }

}

}