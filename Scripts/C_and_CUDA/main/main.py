import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt

def executeForward(x, y, num_segments = 1, tau = 1, ell = 1, p=1, covfunc = "matern", protein_structure = "demoleus2x2"): # "Retinin2x2" or "demoleus2x2"
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

    nx = len(x);
    ny = len(y);

    Xmesh, Ymesh = np.meshgrid(x,y)
    nx, ny = Xmesh.shape
    #print(nx,ny)
    #print(Xmesh)
    #print(Ymesh)
    Xmesh = Xmesh.flatten()
    Ymesh = Ymesh.flatten()

    # Prepare observation points
    n = len(Xmesh)
    x_arr = (ctypes.c_double * n)(*Xmesh)
    y_arr = (ctypes.c_double * n)(*Ymesh)
    

    so_file = "./forward.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    c_func.executeForward(x_arr, y_arr, n, protein_structure_encoded, num_segments)

    # Get result from 
    filename = '../../../Results/forward/Ez_scat.txt'
    data = np.loadtxt(filename)
    Ez_scat_real = data[:,0];
    Ez_scat_imag = data[:,1];
    Ez_scat_real = np.reshape(Ez_scat_real,[nx,ny])
    Ez_scat_imag = np.reshape(Ez_scat_imag,[nx,ny])
    #print(data)

    #print(Ez_scat_real)
    plt.figure()
    plt.imshow(np.sqrt(Ez_scat_real**2 + Ez_scat_imag**2))
    plt.colorbar()
    plt.savefig('Ez_scat.png')

    # Get result from 
    filename = '../../../Results/forward/Hx_scat.txt'
    data = np.loadtxt(filename)
    real = data[:,0];
    imag = data[:,1];
    real = np.reshape(real,[nx,ny])
    imag = np.reshape(imag,[nx,ny])
    #print(data)

    #print(real)
    plt.figure()
    plt.imshow(np.sqrt(real**2 + imag**2))
    plt.colorbar()
    plt.savefig('Hx_scat.png')

    # Get result from 
    filename = '../../../Results/forward/Hy_scat.txt'
    data = np.loadtxt(filename)
    real = data[:,0];
    imag = data[:,1];
    real = np.reshape(real,[nx,ny])
    imag = np.reshape(imag,[nx,ny])
    #print(data)

    #print(real)
    plt.figure()
    plt.imshow(np.sqrt(real**2 + imag**2))
    plt.colorbar()
    plt.savefig('Hy_scat.png')

    # Get result from 
    filename = '../../../Results/forward/Ez_inc.txt'
    data = np.loadtxt(filename)
    real = data[:,0];
    imag = data[:,1];
    real = np.reshape(real,[nx,ny])
    imag = np.reshape(imag,[nx,ny])
    #print(data)

    #print(real)
    plt.figure()
    plt.imshow(np.sqrt((Ez_scat_real + real)**2 + (Ez_scat_imag+imag)**2))
    plt.colorbar()
    plt.savefig('Ez_inc.png')


    plt.figure()
    for i in range(num_segments):
        filename = f'../../../Data/segments/test_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1])
    plt.savefig('test_points.png')

    plt.figure()
    for i in range(num_segments):
        filename = f'../../../Data/segments/ext_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1])
    plt.savefig('ext_points.png')

    plt.figure()
    for i in range(num_segments):
        filename = f'../../../Data/segments/int_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1])
    plt.savefig('int_points.png')

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
    obs_grid = 300;
    Y = np.linspace(0,21*10**(-7),obs_grid);
    Y = Y + 4.38059442329516e-08;
    #Y = Y + 10**(-2);
    #Y = np.linspace(0,10**(-7),obs_grid);
    X = np.linspace(-0.5*10**(-7),20.5*10**(-7),obs_grid);
    
    

    Z = executeForward(x = X, y = Y, num_segments = 10)
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
    
    

    
    


