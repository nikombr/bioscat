
import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import sys

def executeProfileForward(total_grid_points=100, obs_grid=300, num_segments = 1, protein_structure = "demoleus2x2", beta = 0, lambd = 325e-9, deviceComputation = False, location = "near",savefolder=''): # "Retinin2x2" or "demoleus2x2"

    x = np.linspace(-25, 25, obs_grid)*10**(-7);
    y = np.linspace(0,   50, obs_grid)*10**(-7);


    if total_grid_points % num_segments != 0: # Probably not necessary as we will probably only look at one segment
        print("Please check grid and number of segments so number of grid points is the same in every segment!")
        return

    nx = len(x);
    ny = len(y);

    Xmeshsq, Ymeshsq = np.meshgrid(x,y)
    nx, ny = Xmeshsq.shape
    Xmesh = Xmeshsq.flatten()
    Ymesh = Ymeshsq.flatten()

    # Prepare observation points
    n = len(Xmesh)
    x_arr = (ctypes.c_double * n)(*Xmesh)
    y_arr = (ctypes.c_double * n)(*Ymesh)
    
    so_file = "./so/profileForward.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    c_func.executeProfileForward(x_arr, y_arr, n, protein_structure_encoded, num_segments,total_grid_points, ctypes.c_double(beta*math.pi/180), ctypes.c_double(lambd), int(deviceComputation))

if __name__ == "__main__":

    if len(sys.argv) > 1:
        dev = int(sys.argv[1])
        if dev == 1:
            deviceComputation = True
        else:
            deviceComputation = False
            num_threads = int(sys.argv[2])
    else:
        deviceComputation = True;

    dir = "../../../../../../../work3/s194146/bioscatdata"
    if deviceComputation:
        filename = f'{dir}/Results/profiling/grid/forward_device.txt'
        if os.path.isfile(filename):
            os.remove(filename)
    else:
        filename = f'{dir}/Results/profiling/grid/forward_host_{num_threads}.txt'
        if os.path.isfile(filename):
            os.remove(filename)
    obs_grid = 300;
    
    for total_grid_points in range(50,1050,50):
        executeProfileForward(total_grid_points=total_grid_points,obs_grid=obs_grid,deviceComputation=deviceComputation)