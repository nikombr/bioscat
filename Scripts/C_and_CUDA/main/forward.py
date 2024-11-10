import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math

def executeForward(x, y, total_grid_points=100,num_segments = 1, protein_structure = "demoleus2x2", beta = 0, lambd = 325e-9): # "Retinin2x2" or "demoleus2x2"


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
    

    so_file = "./so/forward.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    c_func.executeForward(x_arr, y_arr, n, protein_structure_encoded, num_segments,total_grid_points, ctypes.c_double(beta*math.pi/180), ctypes.c_double(lambd))

    variables = ['x', 'y', 'z'];
    types = ['scat', 'inc', 'ref']
    field_types = ['E','H']

    mdic = dict();

    for j in range(3):
        typ = types[j]
        for k in range(2):
            field_typ = field_types[k]
            fields = np.zeros((nx,ny,3),dtype=complex)

            for i in range(3):
                var = variables[i]
                filename = f'../../../Results/forward/{field_typ}{var}_{typ}.txt'
                data = np.loadtxt(filename)
                field = data[:,0] + 1j*data[:,1]
                fields[:,:,i] = np.reshape(field,[nx,ny])

                plt.figure()
                plt.imshow(np.abs(fields[:,:,i]))
                plt.colorbar()
                plt.savefig(f'plots/{field_typ}{var}_{typ}.png')
                plt.close()
                #time.sleep(30);
                #os.remove(filename)

            mdic[f'{field_typ}_{typ}'] = fields

    savename = f'../../../Results/forward/2D/{protein_structure}/fields_{int(beta)}_lambda_{int(lambd*10**9)}_num_segments_{int(num_segments)}_total_grid_points_{int(total_grid_points)}.mat'
    savemat(savename, mdic)


  

    plt.figure()
    for i in range(num_segments):
        filename = f'../../../Data/segments/test_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1])
    plt.savefig('plots/test_points.png')

    plt.figure()
    for i in range(num_segments):
        filename = f'../../../Data/segments/ext_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1])
    plt.savefig('plots/ext_points.png')

    plt.figure()
    for i in range(num_segments):
        filename = f'../../../Data/segments/int_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        plt.plot(data[:,0],data[:,1])
    plt.savefig('plots/int_points.png')

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
    obs_grid = 200;
    Y = np.linspace(0,21*10**(-7),obs_grid);
    Y = Y + 4.38059442329516e-08;
    #Y = Y + 10**(-2);
    #Y = np.linspace(0,10**(-7),obs_grid);
    X = np.linspace(-0.5*10**(-7),20.5*10**(-7),obs_grid);
    
    

    Z = executeForward(x = X, y = Y, num_segments = 4,beta = 30,total_grid_points=1000)
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
    
    

    
    


