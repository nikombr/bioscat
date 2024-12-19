import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import sys

def executeForward(ystart=-10,total_grid_points=100,num_segments = 1, protein_structure = "demoleus2x2", beta = 0, lambd = 325e-9, deviceComputation = False, location = "near",obs_grid=300,savefolder='',printOutput=True): # "Retinin2x2" or "demoleus2x2"

    x = np.linspace(-25,25,obs_grid)*10**(-7);
    if location == "near":
        y = np.linspace(ystart,ystart+50,obs_grid)*10**(-7);
    elif location == "far":
        y = np.linspace(0,50,obs_grid)*10**(-7);
        y += 3e-2;
    elif location == "far_field_pattern":
        phi = np.linspace(0,np.pi,obs_grid)
        x = 2e-3*np.cos(phi)
        y = 2e-3*np.sin(phi)


    if total_grid_points % num_segments != 0:
        print("Please check grid and number of segments so number of grid points is the same in every segment!")
        return

    nx = len(x);
    ny = len(y);
    if location == "far_field_pattern":
        Xmesh = x;
        Ymesh = y;
    else:
        Xmeshsq, Ymeshsq = np.meshgrid(x,y)
        nx, ny = Xmeshsq.shape
        #print(nx,ny)
        #print(Xmesh)
        #print(Ymesh)
        Xmesh = Xmeshsq.flatten()
        Ymesh = Ymeshsq.flatten()

    # Prepare observation points
    n = len(Xmesh)
    x_arr = (ctypes.c_double * n)(*Xmesh)
    y_arr = (ctypes.c_double * n)(*Ymesh)
    

    so_file = "./so/forward.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    c_func.executeForward(x_arr, y_arr, n, protein_structure_encoded, num_segments,total_grid_points, ctypes.c_double(beta*math.pi/180), ctypes.c_double(lambd), int(deviceComputation), int(printOutput))

    variables = ['x', 'y', 'z'];
    types = ['scat', 'inc', 'int']
    field_types = ['E','H']

    mdic = dict();
    mdic['x'] = x
    mdic['y'] = y
    dir = "../../../../../../../work3/s194146/bioscatdata"

    if not location == "far_field_pattern":
        nanox = f"{dir}/Data/nanostructures/2D/{protein_structure}_x_{total_grid_points}.txt"
        nanox = np.loadtxt(nanox,delimiter=',')
        nanof = f"{dir}/Data/nanostructures/2D/{protein_structure}_f_{total_grid_points}.txt"
        nanof = np.loadtxt(nanof,delimiter=',')
        f = np.interp(x, nanox, nanof)
        f[x > max(nanox)] = 0
        f[x < min(nanox)] = 0
        plt.figure()
        plt.plot(x,f)
        plt.savefig(f'{dir}/tmpplots/interior_f.png',dpi=300)
        plt.close()
        interior = np.zeros((nx,ny))
        for j in range(nx):
            for i in range(ny):
                if Ymeshsq[j,i] <= f[i] and Ymeshsq[j,i] >= 0:
                    interior[j,i] = 1

        plt.figure()
        plt.imshow(interior)
        plt.savefig(f'{dir}/tmpplots/interior.png',dpi=300)
        plt.close()

        bool_compute_interior = ystart < np.max(f) and location == 'near'
        print('int' if bool_compute_interior else 'no')
    else:
        bool_compute_interior = False
    

    for j in range(3 if bool_compute_interior else 2):
        typ = types[j]
        for k in range(2):
            field_typ = field_types[k]
            if location == "far_field_pattern":
                fields = np.zeros((nx,3),dtype=complex)
            else:
                fields = np.zeros((nx,ny,3),dtype=complex)

            for i in range(3):
                var = variables[i]
                filename = f'{dir}/Results/forward/{field_typ}{var}_{typ}.txt'
                data = np.loadtxt(filename)
                field = data[:,0] + 1j*data[:,1]
            
                if location == "far_field_pattern":
                    fields[:,i] = field

                    if typ == 'scat' and var == "z" and field_typ == 'E':
                        plt.figure()
                        print(field.shape)
                        plt.polar(phi,np.abs(field))
                        plt.savefig(f'{dir}/tmpplots/farfieldtest.png',dpi=300)
                        plt.close()
                else:
                    fields[:,:,i] = np.reshape(field,[nx,ny])

                    if typ == 'scat' or typ == 'inc':
                        fields[interior == 1] = 0
                    elif typ == 'int':
                        fields[interior == 0] = 0


                    plt.figure()
                    plt.imshow(np.abs(fields[:,:,i]))
                    plt.colorbar()
                    plt.savefig(f'{dir}/tmpplots/{field_typ}{var}_{typ}.png',dpi=300)
                    plt.close()
                    os.remove(filename)
                    
          
            mdic[f'{field_typ}_{typ}'] = fields
    if location == 'far':
        mdic[f'E_int'] = fields*0
        mdic[f'H_int'] = fields*0
    savename = f'{dir}/Results/forward{savefolder}/{protein_structure}/{location}/fields_beta_{int(beta)}_lambda_{int(lambd*10**9)}_num_segments_{int(num_segments)}_total_grid_points_{int(total_grid_points)}_obs_grid_{int(obs_grid)}.mat'
    savemat(savename, mdic)

    n = 50;
    plt.figure()
    if not location == "far_field_pattern":
        nanox = f"{dir}/Data/nanostructures/2D/{protein_structure}_x_1000.txt"
        nanox = np.loadtxt(nanox,delimiter=',')
        nanof = f"{dir}/Data/nanostructures/2D/{protein_structure}_f_1000.txt"
        nanof = np.loadtxt(nanof,delimiter=',')
        plt.plot(nanox,nanof)
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
        #plt.plot(data[:,0],data[:,1],'b.')
    #plt.savefig('plots/ext_points.png')

    #plt.figure()
    for i in range(num_segments):
        filename = f'{dir}/Data/segments/int_segment_{i+1}.txt'
        data = np.loadtxt(filename)
        #plt.plot(data[:,0],data[:,1],'r.')
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    plt.savefig(f'{dir}/tmpplots/all_points.png',dpi=300)
    plt.close()

def segmentTest(num_segments = 1,beta=0,protein_structure='demoleus2x2'):
    obs_grid = 300
    grid_size = 300
    location = "near"; # near, far
    executeForward(ystart = 16,num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True, location = location,obs_grid=obs_grid,savefolder='/segment_test')
    location = "far"; # near, far
    executeForward(num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True, location = location,obs_grid=obs_grid,savefolder='/segment_test',printOutput=False)

def comsolTest(beta=0,protein_structure='demoleus2x2',obs_grid=300,grid_size=1000):
    num_segments = 1
    location = "near"; # near, far
    executeForward(num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True, location = location, obs_grid=obs_grid,savefolder='/comsol_test',printOutput=False)

def farFieldPatternlTest(beta=0,protein_structure='demoleus2x2',obs_grid=300,grid_size=1000):
    num_segments = 1
    location = "far_field_pattern"; # near, far
    executeForward(num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True, location = location, obs_grid=obs_grid,savefolder='/farFieldPattern',printOutput=False)

def gridPointTest(beta=0,protein_structure='demoleus2x2',obs_grid=300,grid_size=1000):
    num_segments = 1
    location = "near"; # near, far
    executeForward(num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True, location = location, obs_grid=obs_grid,savefolder='/grid_point_test',printOutput=False)

if __name__ == "__main__":
    
    #obs_grid = 200;
    #Y = np.linspace(-10*10**(-7),31*10**(-7),obs_grid);
    #location = "near"; # near, far, (far_field_pattern)
    #Y = Y + 1.6e-6;
    #Y = Y + 3e-2;
    #X = np.linspace(-20.5*10**(-7),20.5*10**(-7),obs_grid);
    #protein_structure = "demoleus2x2"  # "Retinin2x2" or "demoleus2x2"
    
    #executeForward(x = X, y = Y, num_segments = 2, beta = 0, total_grid_points=300, protein_structure = protein_structure, deviceComputation = True, location = location)
    if len(sys.argv) > 1:
        typeTest = int(sys.argv[1])
        
    else:
        typeTest = 0;
        print("Running with default.")
        num_segments = 1
        beta = 0
        obs_grid = 300;
        grid_size=10
        protein_structure = "sintest" # "Retinin2x2" or "demoleus2x2"
        location = "near"
        executeForward(num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True, location = location, obs_grid=obs_grid,savefolder='')
    if typeTest == 1:
        num_segments = int(sys.argv[2])
        beta = int(sys.argv[3])
        protein_structure = sys.argv[4] # "Retinin2x2" or "demoleus2x2"
        segmentTest(num_segments=num_segments, beta = beta, protein_structure = protein_structure)
    elif typeTest == 2:
        beta = int(sys.argv[2])
        protein_structure = sys.argv[3] # "Retinin2x2" or "demoleus2x2"
        obs_grid = int(sys.argv[4])
        comsolTest(beta=beta,protein_structure=protein_structure,obs_grid=obs_grid)
    elif typeTest == 3:
        beta = int(sys.argv[2])
        protein_structure = sys.argv[3] # "Retinin2x2" or "demoleus2x2"
        farFieldPatternlTest(beta=beta,protein_structure=protein_structure)
    elif typeTest == 4:
        beta = int(sys.argv[2])
        protein_structure = sys.argv[3] # "Retinin2x2" or "demoleus2x2"
        grid_size = int(sys.argv[4])
        gridPointTest(beta=beta,protein_structure=protein_structure,grid_size=grid_size)


