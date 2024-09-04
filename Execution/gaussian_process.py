import ctypes
import numpy as np


def gaussian_process(x, y, so_file, tau = 1, ell = 1):
    
    # Compute mesh
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    
    # Get variables for C implementation
    n = len(X)
    X_arr = (ctypes.c_int * len(X))(*X)
    Y_arr = (ctypes.c_int * len(Y))(*Y)

    # Get C function
    c_func = ctypes.CDLL(so_file)

    # Execute C implementation
    c_func.gaussian_process(X_arr, Y_arr, n, tau, ell)

    # Get result from file
    z = 0

    # Return reshaped result
    return np.reshape(z,(len(x),len(y)))