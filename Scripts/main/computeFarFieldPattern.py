import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import sys

def executeFarFieldPattern(total_grid_points=100,num_segments = 1, protein_structure = "demoleus2x2", beta = 0, lambd = 325e-9, deviceComputation = False, savefolder='.',printOutput=True): # "Retinin2x2" or "demoleus2x2"
    if savefolder != "data_test":
        obs_grid = 500;
        phi = np.linspace(0,np.pi,obs_grid)
    else:
        phi = np.array([np.pi/2])

    # Prepare observation points
    n = len(phi)
    phi_arr = (ctypes.c_double * n)(*phi)
    
    so_file = "./so/computeFarFieldPattern.so"

    # Get C function
    c_func = ctypes.CDLL(so_file)
    protein_structure_encoded = protein_structure.encode('utf-8')

    # Execute C implementation
    if savefolder != "data_test":
        c_func.computeFarFieldPattern(phi_arr, n, protein_structure_encoded, num_segments,total_grid_points, ctypes.c_double(beta*math.pi/180), ctypes.c_double(lambd), int(deviceComputation), int(printOutput))

        mdic = dict();
        mdic['phi'] = phi

        dir = "../../../../../../../work3/s194146/bioscatdata"

        filename = f'{dir}/Results/forward/farFieldPattern.txt'
        far_field_pattern = np.loadtxt(filename)
        print(far_field_pattern.shape)
        plt.figure()
        plt.polar(phi, far_field_pattern)
        plt.savefig(f'{dir}/tmpplots/farFieldPattern.png')
        plt.close()
        #os.remove(filename)

        mdic['far_field_pattern'] = far_field_pattern

        location = 'far_field_pattern';
        savename = f'{dir}/Results/forward/{savefolder}/{protein_structure}/{location}/absolute_field_beta_{int(beta)}_lambda_{int(lambd*10**9)}_num_segments_{int(num_segments)}_total_grid_points_{int(total_grid_points)}.mat'
        savemat(savename, mdic)


        plt.figure()
        for i in range(num_segments):
            filename = f'{dir}/Data/segments/test_segment_{i+1}.txt'
            data = np.loadtxt(filename)
            plt.plot(data[:,0],data[:,1],'k.')
        #plt.savefig('plots/test_points.png')

        #plt.figure()
        for i in range(num_segments):
            filename = f'{dir}/Data/segments/ext_segment_{i+1}.txt'
            data = np.loadtxt(filename)
            plt.plot(data[:,0],data[:,1],'.')
        #plt.savefig('plots/ext_points.png')

        #plt.figure()
        for i in range(num_segments):
            filename = f'{dir}/Data/segments/int_segment_{i+1}.txt'
            data = np.loadtxt(filename)
            plt.plot(data[:,0],data[:,1],'.')
        plt.savefig(f'{dir}/tmpplots/all_points.png')
    else:
        dir = "../../../../../../../work3/s194146/bioscatdata"
        filename = f'{dir}/Results/forward/farFieldPattern.txt'
        lambdas = np.linspace(250,750,2000)*1e-9;
        betas = np.arange(0,100,10)
        reflectance = np.zeros((len(lambdas),len(betas)))
        for i, lambd in enumerate(lambdas):
            for j, beta in enumerate(betas):
                print(j+i*len(betas)+1,"/",len(lambdas)*len(betas),flush=True)
                c_func.computeFarFieldPattern(phi_arr, n, protein_structure_encoded, num_segments,total_grid_points, ctypes.c_double(beta*math.pi/180), ctypes.c_double(lambd), int(deviceComputation), int(printOutput))
                reflectance[i,j] = np.loadtxt(filename)
                

        os.remove(filename)
        mdic = dict();
        mdic['lambdas'] = lambdas
        mdic['betas'] = betas
        mdic['reflectance'] = reflectance

        savename = f'{dir}/Results/forward/{savefolder}/{protein_structure}.mat'
        savemat(savename, mdic)





def segmentTest(num_segments = 1,beta=0,protein_structure='demoleus2x2'):
    grid_size = 300
    executeFarFieldPattern(num_segments = num_segments, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True,savefolder='segment_test',printOutput=False)


def comsolTest(beta=0,protein_structure='demoleus2x2',grid_size=1000,lambd=325):
    executeFarFieldPattern(num_segments = 1, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True,savefolder='comsol_test',lambd=lambd,printOutput=False)

def gridPointTest(beta=0,protein_structure='demoleus2x2',grid_size=1000,lambd=325):
    executeFarFieldPattern(num_segments = 1, beta = beta, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True,savefolder='grid_point_test',lambd=lambd,printOutput=False)

def dataTest(protein_structure='demoleus2x2',grid_size=1000):
    executeFarFieldPattern(num_segments = 1, total_grid_points=grid_size, protein_structure = protein_structure, deviceComputation = True,savefolder='data_test',printOutput=False)


if __name__ == "__main__":
    
    
    #location = "far_field_pattern"
    
    #protein_structure = ["demoleus2x2", "Retinin2x2"] # "Retinin2x2" or "demoleus2x2"

    #executeFarFieldPattern(phi=phi, num_segments = 1, beta = 0, total_grid_points=300, protein_structure = protein_structure[0], deviceComputation = True)
    if len(sys.argv) > 1:
        typeTest = int(sys.argv[1])
        
    else:
        print("Please choose something")
    if typeTest == 1:
        num_segments = int(sys.argv[2])
        beta = int(sys.argv[3])
        protein_structure = sys.argv[4] # "Retinin2x2" or "demoleus2x2"
        segmentTest(num_segments=num_segments, beta = beta, protein_structure = protein_structure)
    elif typeTest == 2:
        beta = int(sys.argv[2])
        protein_structure = sys.argv[3] # "Retinin2x2" or "demoleus2x2"
        lambd = int(sys.argv[4])*1e-9
        print("lambda = ",lambd)
        comsolTest(beta=beta,protein_structure=protein_structure,lambd=lambd)
    elif typeTest == 4:
        beta = int(sys.argv[2])
        protein_structure = sys.argv[3] # "Retinin2x2" or "demoleus2x2"
        lambd = int(sys.argv[4])*1e-9
        grid_size = int(sys.argv[5])
        gridPointTest(beta=beta,protein_structure=protein_structure,lambd=lambd,grid_size=grid_size)
    elif typeTest == 5:
        protein_structure = sys.argv[2] # "Retinin2x2" or "demoleus2x2"
        dataTest(protein_structure=protein_structure)



