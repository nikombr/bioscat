#include <stdlib.h>
#include <stdio.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string.h>
extern "C" {
#include "../../lib/Segment.h"
#include "../../lib/RealMatrix.h"
using namespace std;

void Segment::computeFieldsForLinearSystem() {

    
    computeIncidentFieldVectors(y_test);
    //printf("hej1\n");
    computeReflectedFieldVectors(y_test);
    //printf("hej2\n");
    computeScatteredFieldMatrices(x_test, y_test, false);
    computeInteriorFieldMatrices(x_test, y_test);

}

}