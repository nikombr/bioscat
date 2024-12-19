import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import shutil
import sys

def executeInverse(num_segments = 1, delta = 0.01, maxiter = 10, total_grid_points=100,protein_structure = "demoleus2x2",tau = 1, ell = 1, p=1, covfunc = "squared_exponential", decay_rate = 0, gamma=1e5, fine_tuning = False,datatype='angle_resolved',testTypeName=None): # "Retinin2x2" or "demoleus2x2"    
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
    if testTypeName == None:
        filename = f"{protein_structure}/{covfunc}/{spec}_delta_{delta}_gamma_{gamma:.0e}".replace("+0", "")
    else:
        filename = f"{testTypeName}/{protein_structure}_{spec}_delta_{delta}_gamma_{gamma:.0e}".replace("+0", "")
    hyper_arr = (ctypes.c_double * len(hyper))(*hyper)

    so_file = "./so/inverse.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    c_func.inverse(protein_structure_encoded, num_segments, total_grid_points, hyper_arr, len(hyper_arr), type_covfunc, ctypes.c_double(delta), maxiter, filename.encode('utf-8'), ctypes.c_double(decay_rate), ctypes.c_double(gamma), int(fine_tuning), datatype.encode('utf-8'))
    dir = "../../../../../../../work3/s194146/bioscatdata"

    plt.figure()
    for i in range(num_segments):
        filename = f'{dir}/Data/segments/test_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1],'k.')
    for i in range(num_segments):
        filename = f'{dir}/Data/segments/n_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        #plt.plot(data[:10,0],data[:10,1],'k.')
    #plt.savefig('plots/test_points.png')

    #plt.figure()
    for i in range(num_segments):
        filename = f'{dir}/Data/segments/ext_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1],'b.')
    #plt.savefig('plots/ext_points.png')

    #plt.figure()
    for i in range(num_segments):
        filename = f'{dir}/Data/segments/int_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1],'r.')
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig(f'{dir}/tmpplots/all_points.png',dpi=300)
    plt.close()

def maternParameterTest(protein_structure, datatype, p = 2, tau = 0.5, ell = 0.07, delta = 0.005, gamma = 1e5):
    executeInverse(testTypeName = 'maternParameterTest', p = p, tau = tau, ell = ell, gamma = gamma, delta = delta, datatype = datatype, total_grid_points=300, maxiter=20000, covfunc = "matern", protein_structure = protein_structure, decay_rate=0, fine_tuning = False) 


def squaredExponentialParameterTest(protein_structure, datatype, tau = 0.5, ell = 0.07, delta = 0.005, gamma = 1e5):
    executeInverse(testTypeName = 'squaredExponentialParameterTest', tau = tau, ell = ell, gamma = gamma, delta = delta, datatype = datatype, total_grid_points=300, maxiter=20000, covfunc = "squared_exponential", protein_structure = protein_structure, decay_rate=0, fine_tuning = False) 


if __name__ == "__main__":
    if len(sys.argv) > 1:
        typeTest = int(sys.argv[1])
    else:
        typeTest = 0
        # Rigtig god indstilling it seems:
        # executeInverse(p=1,ell=0.5,total_grid_points=300,num_segments=1,delta=0.1,maxiter=20000,covfunc = "matern",protein_structure = "Retinin2x2",decay_rate=1e-3, fine_tuning = True, gamma=1e4) # squared_exponential matern
        #executeInverse(p=1,ell=0.5,total_grid_points=300,num_segments=1,delta=0.1,maxiter=20000,covfunc = "matern",protein_structure = "Retinin2x2",decay_rate=0.1, fine_tuning = True, gamma=1e4) # squared_exponential matern
        executeInverse(p=4,tau =0.5,ell=0.07,total_grid_points=200,delta=0.01,maxiter=20000,covfunc = "matern",protein_structure = "Retinin2x2",decay_rate=0, fine_tuning = False, gamma=5e4, datatype='angle_resolved') # squared_exponential matern
        #executeInverse(tau =0.7,ell=0.07,total_grid_points=200,delta=0.001,maxiter=20000,covfunc = "squared_exponential",protein_structure = "demoleus2x2",decay_rate=0, fine_tuning = False, gamma=1e4, datatype='angle_resolved')

    if typeTest == 1: # matern, testing different parameter settings, no tuning or decay, different p, tau, ell
        protein_structure = sys.argv[2] # "Retinin2x2" or "demoleus2x2"
        datatype          = sys.argv[3] # 'angle_resolved' or 'one_observation'
        p                 = int(sys.argv[4])
        tau               = float(sys.argv[5])
        ell               = float(sys.argv[6])
        maternParameterTest(protein_structure, datatype, p = p, tau = tau, ell = ell)

    elif typeTest == 2: # squared exponential, testing different settings, no tuning or decay, different tau, ell, delta, gamma
        protein_structure = sys.argv[2] # "Retinin2x2" or "demoleus2x2"
        datatype          = sys.argv[3] # 'angle_resolved' or 'one_observation'
        tau               = float(sys.argv[4])
        ell               = float(sys.argv[5])
        squaredExponentialParameterTest(protein_structure, datatype, tau = tau, ell = ell)
    elif typeTest == 3: # matern, use best settings for parameters and a bit surrounding perhaps, try different decay rate, try different gamma
        pass
    elif typeTest == 4: # squared exponential, use best settings for parameters and a bit surrounding perhaps, try different decay rate, try different gamma
        pass
    elif typeTest == 5: # matern, combine everything
        pass
    elif typeTest == 6: # squared exponential, combine everything
        pass
    
    
    

    
    


