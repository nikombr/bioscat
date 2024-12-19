import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import shutil

datatype = 'angle_resolved' # 'one_observation'  'angle_resolved'
protein_structure = 'Retinin2x2';
total_grid_points = 1000
testTypeName = 'maternParameterTest' # maternParameterTest squaredExponentialParameterTest


covfunc = 'matern'; # squared_exponential matern
p = 2
tau = 0.7
ell = 0.07
delta = 0.005
gamma=1e5

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
    

if testTypeName == None:
    filenamespec = f"{protein_structure}/{covfunc}/{spec}_delta_{delta}_gamma_{gamma:.0e}".replace("+0", "")
else:
    filenamespec = f"{testTypeName}/{protein_structure}_{spec}_delta_{delta}_gamma_{gamma:.0e}".replace("+0", "")

dir = "../../../../../../../work3/s194146/bioscatdata"
plt.figure()
filename = f'{dir}/Data/segments/test_segment_1.txt'
data1 = np.loadtxt(filename)

nanox = f"{dir}/Data/nanostructures/2D/{protein_structure}_x_{total_grid_points}.txt"
nanox = np.loadtxt(nanox,delimiter=',')
nanof = f"{dir}/Data/nanostructures/2D/{protein_structure}_f_{total_grid_points}.txt"
nanof = np.loadtxt(nanof,delimiter=',')

#filename = f"{dir}/Results/inverse/output.txt"
if covfunc == "squared_exponential":
    type_covfunc = 1
    hyper = np.array([tau,ell])
    spec = f"tau_{tau}_ell_{ell:.2f}"
elif covfunc == "matern":
    type_covfunc = 2
    hyper = np.array([tau,ell,p])
    spec = f"tau_{tau}_ell_{ell:.2f}_p_{p}"
filename = f"{dir}/Results/inverse/{datatype}/{filenamespec}_output.txt"
data = np.loadtxt(filename);
[m, n] = data.shape
print(f"{m} accepted")
x = np.linspace(np.min(data1[:,0]),np.max(data1[:,0]),n)
for i in range(m):
    plt.plot(x,data[i,:],'m-',linewidth=1, alpha=i/m)
half = int(m/2)
plt.plot(nanox,nanof,'k-')
y = np.mean(data[half:,:],axis=0)
print("mu = ",np.mean(y))
print("std = ",np.std(y))
plt.plot(x,y,'b-')

plt.plot(x,y*0+np.mean(y),'g--')
plt.plot(x,y*0+np.mean(y)+np.std(y),'g--')
plt.plot(x,y*0+np.mean(y)-np.std(y),'g--')


plt.savefig(f"{dir}/tmpplots/inverse_last_temp.png")
plt.close()