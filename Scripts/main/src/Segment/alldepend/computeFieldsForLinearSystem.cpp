#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
#include "Segment.h"
extern "C" {
using namespace std;

void Segment::computeFieldsForLinearSystem() {

    
    computeIncidentFields(test_points.y);
    computeScatteredFieldMatrices(test_points.x, test_points.y);
    computeInteriorFieldMatrices(test_points.x, test_points.y);

}

}