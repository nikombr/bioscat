import ctypes
import numpy as np


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
    z = 0

    # Return reshaped result
    #return np.reshape(z,(len(x),len(y)))
    return z

if __name__ == "__main__":
    x = np.linspace(-1,1,3)
    y = np.linspace(-1,1,3)
    #y = np.array([]);
    Z = gaussian_process(x, y, tau = 1, ell = 1)
    print("Hej fra Python")