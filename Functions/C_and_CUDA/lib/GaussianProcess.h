#ifndef _INVERSE_PLANE_H
#define _INVERSE_PLANE_H
// x_h, y_h: x and y coordinates for measurement points on host
    // n: number of measurement points
    // hyper: array of hyperparameters
    // num: number of hyperparameters
class GaussianProcess {
    private:        
        bool device;        // true if a device is available else false
        int n;              // number of points for estimating plane
        double *hyper_h;    // hyperparameters on host
        double *hyper_d;    // hyperparameters on device
        int num;            // number of hyperparameters
        int dim;            // dimension of the problem, either 1 for curve or 2 for plane
        double *x_h;        // x coordinates for estimating plane on host
        double *x_d;        // x coordinates for estimating plane on device
        double *y_h;        // y coordinates for estimating plane on host
        double *y_d;        // y coordinates for estimating plane on device
        double *z_h;        // height of plane in location (x,y) on host
        double *z_d;        // height of plane in location (x,y) on device
        double *p_h;        // random vector on host
        double *p_d;        // random vector on device
        double **Sigma_h;   // covariance matrix on host
        double **Sigma_d;   // covariance matrix on device
        double *Sigma_log;  // Sigma_d[0] on device
        double **L_h;       // lower triangular matrix from Cholesky factorization on host
        double **L_d;       // lower triangular matrix from Cholesky factorization on device
        double *L_log;      // L_d[0] on device

    public:
        GaussianProcess(double* x_h, double* y_h, int n, double* hyper, int num, int dim);  // Constructer, sets default values and allocates
        ~GaussianProcess();                                                                 // Destructer
        void covariance_matrix();                                                           // Computes covariance matrix K
        void cholesky();                                                                    // Does cholesky factorization of K to compute L
        void realisation();                                                                 // Computes realisation of the Gaussian process from L
};

#endif