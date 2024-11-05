import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import shutil

def executeGenerateArtificialData(num_obs_points=30, num_segments = 1, total_grid_points=10,protein_structure = "demoleus2x2"): # "Retinin2x2" or "demoleus2x2"

    phi = np.linspace(0,math.pi,num_obs_points);
    x = 10**(-2)*np.cos(phi)
    y = 10**(-2)*np.sin(phi)

    directory_name = f'../../../Data/artificial_data/{protein_structure}/run_num_obs_points_{num_obs_points}_num_segments_{num_segments}_total_grid_points_{total_grid_points}/'
    if os.path.exists(directory_name):
        shutil.rmtree(directory_name)
        os.makedirs(directory_name)  # Recreate the directory if you want to keep it
    else:
        os.mkdir(directory_name)

    # Prepare observation points
    n = len(x)
    x_arr = (ctypes.c_double * n)(*x)
    y_arr = (ctypes.c_double * n)(*y)
    

    so_file = "./so/generateArtificialData.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    lambda0 = 325e-9;
    lambdas = np.linspace(0.5*lambda0,1.5*lambda0,3);
    betas = np.linspace(0,math.pi,3);

    lambdas_arr = (ctypes.c_double * len(lambdas))(*lambdas)
    betas_arr = (ctypes.c_double * len(betas))(*betas)

    # Execute C implementation
    c_func.executeGenerateArtificialData(x_arr, y_arr, n, protein_structure_encoded, num_segments, total_grid_points,betas_arr, lambdas_arr,len(betas),len(lambdas))

    # Find files and move them to the correct directory
    source = '../../../Data/artificial_data/temp/reflectance.txt'
    shutil.move(source, directory_name)
    
    np.savetxt(f'{directory_name}x_obs.txt', x, fmt='%e', delimiter='\n')
    np.savetxt(f'{directory_name}y_obs.txt', y, fmt='%e', delimiter='\n')
    np.savetxt(f'{directory_name}lambdas.txt', lambdas, fmt='%e', delimiter='\n')
    np.savetxt(f'{directory_name}betas.txt', betas, fmt='%e', delimiter='\n')



    data = 0;
    return data

if __name__ == "__main__":
    
    

    Z = executeGenerateArtificialData()
    """
    A_real = np.loadtxt("A_real.txt")
    A_real_C = np.loadtxt("A_real_C.txt")
    diff = A_real-A_real_C
    plt.figure()
    plt.imshow(np.log(np.abs(diff)/np.abs(A_real)))
    plt.colorbar()
    plt.savefig('real.png')
    print(f"error real A = {np.max(diff)}")

    A_imag = np.loadtxt("A_imag.txt")
    A_imag_C = np.loadtxt("A_imag_C.txt")
    diff = A_imag-A_imag_C
    plt.figure()
    plt.imshow(np.abs(diff)/np.abs(A_imag))
    plt.colorbar()
    plt.savefig('imag.png')
    print(f"error imag A = {np.max(diff)}")


    b_real = np.loadtxt("b_real.txt")
    b_real_C = np.loadtxt("b_real_C.txt")
    diff = b_real-b_real_C
    print(f"error real b = {np.max(diff)}")

    b_imag = np.loadtxt("b_imag.txt")
    b_imag_C = np.loadtxt("b_imag_C.txt")
    diff = b_imag-b_imag_C
    print(f"error imag b = {np.max(diff)}")
    

    bbig = np.loadtxt("bbig.txt")
    bbig_C = np.loadtxt("bbig_C.txt")
    Abig = np.loadtxt("Abig.txt")
    Abig_C = np.loadtxt("Abig_C.txt")
    diff = np.abs(bbig-bbig_C)
    print(f"error bbig = {np.max(diff)}")
    diff = np.abs(Abig-Abig_C)
    print(f"error Abig = {np.max(diff)}")
    """
    
    

    
    


