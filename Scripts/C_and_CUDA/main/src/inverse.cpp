#include "../lib/forward.h"

void executeInverse(double *x, double*y, int n ,char* protein_structure, int num_segments, double beta, double lambda) {

    forward(x, y, n, protein_structure, num_segments);

}