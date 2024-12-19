import ctypes
import numpy as np
from scipy.io import savemat
import os
import time

def gaussian_process(x, y, tau = 1, ell = 1, p=1, device = True, covfunc = "squared_exponential"):
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
    

    so_file = "./so/computeGaussianProcess.so"
    
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
    c_func.computeGaussianProcess(X_arr, Y_arr, n, hyper_arr, len(hyper_arr), dim, dev, type_covfunc)

    # Get result from file
    dir = "../../../../../../../work3/s194146/bioscatdata"
    filename = f'{dir}/Results/gaussian_process_realisations/output.txt'
    data = np.loadtxt(filename)
    #print(data)

    # Remove txt output
    os.remove(filename)

    # Save in mat file
    if len(y) == 0:
        mdic = {"data": data, "x": x}
        savemat(f"{dir}/Results/gaussian_process_realisations/curve/{covfunc}/realisations_{spec}.mat", mdic)
    else:
        mdic = {"data": data, "x": x, "y": y}
        savemat(f"{dir}/Results/gaussian_process_realisations/plane/{covfunc}/realisations_{spec}.mat", mdic)
            
    return data

if __name__ == "__main__":
    x = np.linspace(-1,1,300)
    #y = np.linspace(0,10,300)
    y = np.array([]);

    covfunc = "matern" # "squared_exponential" "matern"

    if covfunc == "squared_exponential":
        taus = np.round(np.arange(0.3,1.3,0.2),1);
        ells = np.round(np.arange(0.04,0.13,0.03),2);
        for tau in taus:
            for ell in ells:
                print(tau, ell)
                Z = gaussian_process(x, y, tau = tau, ell = ell, device = True, covfunc = covfunc)
                time.sleep(0.1)

    elif covfunc == "matern":

        taus = np.round(np.arange(0.3,1.3,0.2),1);
        ells = np.round(np.arange(0.04,0.13,0.03),2);
        ps = np.arange(1,5,1);

        for tau in taus:
            for ell in ells:
                for p in ps:
                    print(tau, ell, p)
                    Z = gaussian_process(x, y, p = p, ell = ell, tau = tau, device = True, covfunc = covfunc)
                    time.sleep(0.1)
