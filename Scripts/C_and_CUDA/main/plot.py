import ctypes
import numpy as np
from scipy.io import savemat
import os
import time
import matplotlib.pyplot as plt
import math
import shutil



plt.figure()
filename = f'../../../Data/segments/test_segment_1.txt'
data1 = np.loadtxt(filename)



filename = "../../../Results/inverse/output.txt"
data = np.loadtxt(filename);
[m, n] = data.shape
x = np.linspace(np.min(data1[:,0]),np.max(data1[:,0]),n)
for i in range(m):
    plt.plot(x,data[i,:],'y-',linewidth=1, alpha=i/m)
half = int(m/2)
plt.plot(data1[:,0],data1[:,1],'k-')
y = np.mean(data[half:,:],axis=0)
plt.plot(x,y,'r-')


plt.savefig("plots/inverse_last_temp.png")
plt.close()