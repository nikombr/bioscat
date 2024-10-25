import ctypes
import numpy as np
from scipy.io import savemat
import os
import time

def executeForward(num_segments = 1, tau = 1, ell = 1, p=1, covfunc = "matern", protein_structure = "demoleus2x2"): # "Retinin2x2" or "demoleus2x2"
    # covfunc: covariance function
    if covfunc == "squared_exponential":
        type_covfunc = 1
        hyper = np.array([tau,ell])
        var = "tau"
    elif covfunc == "matern":
        type_covfunc = 2
        hyper = np.array([p,ell])
        var = "p"
        tau = p
    else:
        print("Please input a valid covariance function")
        return
    

    so_file = "./forward.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    c_func.executeForward(protein_structure_encoded, num_segments)

    # Get result from file
    #data = np.loadtxt('../../Data/gaussian_process_realisations/output.txt')
    #print(data)

    # Remove txt output
    """os.remove("../../Data/gaussian_process_realisations/output.txt")

    # Save in mat file
    if len(y) == 0:
        mdic = {"data": data, "x": x}
        savemat(f"../../Data/gaussian_process_realisations/curve_{covfunc}_{var}_{tau}_ell_{ell}.mat", mdic)
    else:
        mdic = {"data": data, "x": x, "y": y}
        savemat(f"../../Data/gaussian_process_realisations/plane_{covfunc}_{var}_{tau}_ell_{ell}.mat", mdic)
"""
    data = 0;
    return data

if __name__ == "__main__":

    Z = executeForward(num_segments = 10)


