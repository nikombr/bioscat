import ctypes
import numpy as np
from scipy.io import savemat
import os

def gaussian_process(x, y, tau = 1, ell = 1):

    so_file = "./gaussian_process.so"
    
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
    hyper = np.array([tau,ell])
    hyper_arr = (ctypes.c_double * len(hyper))(*hyper)

    # Get C function
    c_func = ctypes.CDLL(so_file)

    # Execute C implementation
    c_func.gaussian_process(X_arr, Y_arr, n, hyper_arr, 2, dim)

    # Get result from file
    data = np.loadtxt('../../Data/gaussian_process_realisations/output.txt')
    #print(data)

    # Remove txt output
    os.remove("../../Data/gaussian_process_realisations/output.txt")

    # Save in mat file
    if len(y) == 0:
        mdic = {"data": data, "x": x}
        savemat(f"../../Data/gaussian_process_realisations/curve_squared_exponential_tau_{tau}_ell_{ell}.mat", mdic)
    else:
        mdic = {"data": data, "x": x, "y": y}
        savemat(f"../../Data/gaussian_process_realisations/plane_squared_exponential_tau_{tau}_ell_{ell}.mat", mdic)
            
    return data

if __name__ == "__main__":
    x = np.linspace(0,10,200)
    y = np.linspace(0,10,200)
    #y = np.array([]);
    ells = [1, 2, 4];
    taus = [0.25, 0.5, 1];
    for ell in ells:
        for tau in taus:
            Z = gaussian_process(x, y, tau = tau, ell = ell)
    print("Hej fra Python")