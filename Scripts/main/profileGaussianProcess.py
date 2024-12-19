import ctypes
import numpy as np
from scipy.io import savemat
import os
import time

def profile_gaussian_process(x, y, tau = 1, ell = 1, p=1, device = True, covfunc = "squared_exponential"):
    # covfunc: covariance function
    if covfunc == "squared_exponential":
        type_covfunc = 1
        hyper = np.array([tau,ell])
        spec = f"tau_{tau}_ell_{ell:.2f}"
    elif covfunc == "matern":
        type_covfunc = 2
        hyper = np.array([tau,ell,p])
        spec = f"tau_{tau}_ell_{ell:.2f}_p_{p}"
    else:
        print("Please input a valid covariance function")
        return
    

    so_file = "./so/profileGaussianProcess.so"
    
    # Compute mesh if we are computing planes
    if len(y) != 0:
        X, Y = np.meshgrid(x, y)
        X = X.flatten()
        Y = Y.flatten()
        dim = 2
    else:
        X = x
        Y = y
        dim = 1
    
    # Get variables for C implementation
    n = len(X)
    X_arr = (ctypes.c_double * len(X))(*X)
    Y_arr = (ctypes.c_double * len(Y))(*Y)
    hyper_arr = (ctypes.c_double * len(hyper))(*hyper)
    dev = 1 if device else 0

    # Get C function
    c_func = ctypes.CDLL(so_file)

    # Execute C implementation
    c_func.profileGaussianProcess(X_arr, Y_arr, n, hyper_arr, len(hyper_arr), dim, dev, type_covfunc)



if __name__ == "__main__":
    x = np.linspace(-1,1,300)
    #y = np.linspace(0,10,300)
    y = np.array([]);

    dir = "../../../../../../../work3/s194146/bioscatdata"
    filename = f'{dir}/Results/profiling/GaussianProcess_matern.txt'
    if os.path.isfile(filename):
        os.remove(filename)
    filename = f'{dir}/Results/profiling/GaussianProcess_squared_exponential.txt'
    if os.path.isfile(filename):
        os.remove(filename)

    covfunc = "matern" # "squared_exponential" "matern"
    for n in range(1000,0,-50):
        x = np.linspace(-1,1,n)
        profile_gaussian_process(x, y, p = 1, ell = 1, tau = 1, device = True, covfunc = covfunc)
