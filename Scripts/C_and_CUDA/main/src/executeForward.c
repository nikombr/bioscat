#include <stdlib.h>
#include <stdio.h>
#include "../lib/GP/gaussian_process_inner.h"
#include "../lib/forward.h"
#include <cuda_runtime_api.h>

void executeForward(const char* protein_structure) {

    forward(protein_structure);

}