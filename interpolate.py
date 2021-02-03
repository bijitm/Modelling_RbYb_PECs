import numpy as np
from scipy import interpolate
abinitioFile = np.loadtxt("./PECs_so_fit.dat")
R = abinitioFile[:,0]
print("Original no. of grids:",R.size)
V = abinitioFile[:,1:]
R_new = np.arange(3,4,0.1)
# R_new = np.append(R_new,np.arange(4,8,0.02))
R_new = np.append(R_new,np.arange(4,10,0.05))
R_new = np.append(R_new,np.arange(10,12,0.1))
R_new = np.append(R_new,np.arange(12,20.5,0.5))
print("New no. of grids:",R_new.size)

Vnew = np.zeros([R_new.size,9])
for i in range(9):
    f = interpolate.interp1d(R,V[:,i])
    Vnew[:,i] = f(R_new)

txt = ""
for i in range(R_new.size):
    txt += "{:.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\
\t{:10.2f}\t{:10.2f}\t{:10.2f}\n".format(R_new[i],*Vnew[i,:])

with open("PECs_so_fit_dense.dat","w") as f:
    f.write(txt)

