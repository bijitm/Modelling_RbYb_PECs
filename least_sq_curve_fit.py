import numpy as np
from numpy import random
from numpy import linalg as LA
from scipy.optimize import least_squares


def damp(R,beta,n):
    '''Tang-Toennies damping function'''
    s = 0
    for m in range(n+1):
        s += (np.outer(R,beta))**m/np.math.factorial(m)
    return 1 - np.exp(-np.outer(R,beta))*s

def model_pecs(R,param):
    '''Spin-free electrostatic potentials modeled as:
    Aexp(-beta1*R) - Bexp(-beta2*R) - Sum_{n=6,8}D_n*C_n/R^n
    D_n = 1 - exp(-beta1*R)Sum_{i=0}^n (beta1*R)^i/i! '''

    param = param.reshape(nparam_pec, npec)
    A = param[0,:]*10*Eh
    B = param[1,:]*10*Eh
    beta1 = param[2,:]/a0
    beta2 = param[3,:]/a0
    c6pi, c6sig = 4067.6*Eh*a0**6, 4661.5*Eh*a0**6
    C6 = np.array([c6sig,c6pi,c6pi,c6sig])
    C8 = np.copy(C6)
    C8 = C8*80*a0**2

    #repulsion term
    V_rep = A*np.exp(-np.outer(R,beta1))
    #attraction term
    V_att = -B*np.exp(-np.outer(R,beta2))
    V_att += -damp(R,beta1,6)*np.outer(1/R**6,C6) - damp(R,beta1,8)*np.outer(1/R**8,C8) 

    V = V_rep + V_att
    return V

def so_operator(R,param):
    '''a_SO model functions'''

    param = param.reshape(nparam_aso, naso)
    a_so = 806.09 
    r0 = param[0,:]
    c = param[1,:]*100
    alpha = param[2,:]

    a = a_so + c*(1-np.tanh(alpha*np.add.outer(R,-r0)))
    return a

def func_pecs(vmat,a_so):
    '''SO coupled PECs obtained from model
    spin-free PECs and a_SO functions'''

    w_so = np.zeros([vmat.shape[0],9])
    for i in range(vmat.shape[0]):
        a1, a2, a3, a4 = a_so[i,:]
        v1, v2, v3, v4 = vmat[i,:]
        u1 = np.diag((v3,v2,v4))
        u2 = np.diag((v3,v2,v4,v1,v3))

        u_so1 = np.array([
                        [a3/3,np.sqrt(2*a2*a3)/3,np.sqrt(2*a3*a4/3)],
                        [np.sqrt(2*a2*a3)/3,2*a2/3,-np.sqrt(a2*a4/3)],
                        [np.sqrt(2*a3*a4/3),-np.sqrt(a2*a4/3),0]
                        ])

        u_so2 = np.array([
                [-a3/3,np.sqrt(2*a2*a3)/3,2*np.sqrt(2*a3*a4)/3,np.sqrt(a1*a3)/3,0],
                [np.sqrt(2*a2*a3)/3,-2*a2/3,-np.sqrt(a2*a4)/3,2*np.sqrt(2*a1*a2)/3,0],
                [2*np.sqrt(2*a3*a4)/3,-np.sqrt(a2*a4)/3,0,0,np.sqrt(2*a3*a4/3)],
                [np.sqrt(a1*a3)/3,2*np.sqrt(2*a1*a2)/3,0,0,-np.sqrt(a1*a3/3)],
                [0,0,np.sqrt(2*a3*a4/3),-np.sqrt(a1*a3/3),-a3]
                        ])

        w1,_ = LA.eigh(u1 + u_so1)
        w2,_ = LA.eigh(u2 + u_so2)

        w = np.append(w1,w2)
        w_so[i,:] = np.append([v3+a3],w)
    
    return w_so



def func_residues(x,*args):
    '''Calculates array of residues (y_fit - y_data) for all R grids
    and all 9 SO-coupled PECs. Takes the R grid and the ab inito PEC 
    data as args. x is the array of all parameters (PECs and a_SO).'''

    R_arr, v_so = args
    param_pecs, param_aso = x[:npec*nparam_pec], x[npec*nparam_pec:]
    v = model_pecs(R_arr,param_pecs)         #calls model spin-free potentials
    a_so = so_operator(R_arr,param_aso)      #calls model a_SO functions
    w_so = func_pecs(v,a_so)               #calculates the 9 SO-coupled PECs

    # assigning weights
    #lowering weights for any R<=3.5 ang
    w_so[R_arr<=3.5,:] = 0.01*w_so[R_arr<=3.5,:]
    v_so[R_arr<=3.5,:] = 0.01*v_so[R_arr<=3.5,:]
    
    v_so = v_so.reshape(R_arr.size * 9)
    w_so = w_so.reshape(R_arr.size * 9)

    print(np.sum((v_so-w_so)**2)/R_arr.size)

    return  v_so - w_so




Eh, a0 = 2.1947463e5, 0.5291772         #Hartree to cm-1, bohr to ang
npec, nparam_pec = 4, 4                 #No. of PECs, No. of parameters for each PEC
naso, nparam_aso = 4, 3                 #No. of a_SO functions, parameters for each func

parameters = np.loadtxt("./fortran_pec_input.dat")
param_pecs = parameters[:4,:]
param_aso = parameters[4:,:]

#Guess parameters for spin-free PECs
#Rows ----> parameters,  Columns ----> PECs
# 1st row --> A (in units of Eh)
# 2nd row --> B (in units of Eh)
# 3rd row --> beta1 (in units of a0^-1)
# 4th row --> beta2 (in units of a0^-1)

param_pecs = param_pecs.reshape(param_pecs.size)

#Guess parameters for a_SO functions
#Rows ----> parameters,  Columns ----> 4 functions
# 1st row --> R0 (in units of ang)
# 2nd row --> c (cm^-1)
# 3rd row --> alpha (in units of ang^-1)

param_aso = param_aso.reshape(param_aso.size)

x0 = np.append(param_pecs,param_aso)
ll = nparam_pec*npec

#Defining an array of uncertainties in the parameters
A_sig, beta_sig = 0.001, 0.0001
# A_sig, beta_sig = 0, 0
c_sig, r0_sig, alpha_sig = 0.5, 0.05, 0.5
x_sig = random.rand(x0.size)
x_sig[0:8] = A_sig * (1-2*x_sig[0:8])
x_sig[8:16] = beta_sig * (1-2*x_sig[8:16])
x_sig[ll:ll+4] = r0_sig * (1-2*x_sig[ll:ll+4])
x_sig[ll+4:ll+8] = c_sig * (1-2*x_sig[ll+4:ll+8])
x_sig[ll+8:ll+12] = alpha_sig * (1-2*x_sig[ll+8:ll+12])
x0 = x0 + x_sig

#defining bounds on the parameters
x_lbounds = np.zeros(x0.size)               #default lower bounds for all
x_ubounds = np.inf * np.ones(x0.size)       #default upper bounds for all
# x_lbounds[ll:ll+4] = 4                      #lower bounds for R0
# x_lbounds[ll+4:ll+7] = -1                   #lower bounds for c
x_lbounds[ll+4:ll+7] = -4.03045            #lower bounds for c
x_lbounds[ll+7] = -0.5                           #lower bounds for c
x_lbounds[ll+8:ll+12] = 0.5                     #lower bounds for alpha
x_ubounds[7] = 1.0e-12                        #upper bound for B_^4Sigma
x_ubounds[15] = 1.0e-12                        #upper bound for beta2_^4Sigma
# x_ubounds[ll:ll+4] = 6                      #upper bounds for R0
# x_ubounds[ll+3] = 1.e-12                   #upper bound for R0 of a4
x_ubounds[ll+7] = 0.5                   #upper bound for c of a4
# x_ubounds[ll+11] = 1.e-12                   #upper bound for alpha of a4

for i in range(x0.size):
    if x0[i] < x_lbounds[i]:
        x0[i] = x_lbounds[i]
    if x0[i] > x_ubounds[i]:
        x0[i] = x_ubounds[i]

#Reading interpolated datafile
abinitioFile = np.loadtxt("./PECs_so_fit_dense.dat")
R_arr = abinitioFile[:,0]
print("No. of R grids for optimization:",R_arr.size,"\n")
vv = abinitioFile[:,1:]
v_so = np.zeros(vv.shape)

#Rerranging PEC columns according to decreasing values of omega
#Columns order ---> 5/2, 3/2, 3/2, 3/2, 1/2, 1/2, 1/2, 1/2, 1/2
v_so[:,0] = vv[:,7]
v_so[:,1] = vv[:,2]
v_so[:,2] = vv[:,5]
v_so[:,3] = vv[:,8]
v_so[:,4] = vv[:,0]
v_so[:,5] = vv[:,1]
v_so[:,6] = vv[:,3]
v_so[:,7] = vv[:,4]
v_so[:,8] = vv[:,6]

print("Running optimization...\n")

args = (R_arr, v_so)
res = least_squares(func_residues,x0,bounds=(x_lbounds,x_ubounds),\
                    ftol=1.e-7,xtol=1.e-8,gtol=1.e-8,args=args,verbose=1)
param_pecs, param_aso = res.x[:npec*nparam_pec], res.x[npec*nparam_pec:]

print("Optimization done\n")

param_wr = np.copy(res.x)

txt1 = "{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n\
{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n\
{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n\
{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n\n\
{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n\
{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n\
{:18.6f}\t{:18.6f}\t{:18.6f}\t{:18.6f}\n".format(*param_wr)
with open("output_param.dat","w") as f:
    f.write(txt1)

txt2 = "cost: " + "{:20.8f}\n\n".format(res.cost/R_arr.size)
with open("cost.dat","a") as f:
    f.write(txt2)

# Defining new grid space as the original abinitio datafile
R_new = np.arange(3,6,0.05)
R_new = np.append(R_new,np.arange(6,8,0.1))
R_new = np.append(R_new,np.arange(8,9,0.25))
R_new = np.append(R_new,np.arange(9,20.5,0.5))
print("No. of R grids to write data:",R_new.size,"\n")

vv = model_pecs(R_new,param_pecs)
a_so = so_operator(R_new,param_aso)
pecs = func_pecs(vv,a_so)

txt3,txt4 = "",""
for i in range(R_new.size):
    txt3 += "{:.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\
\t{:10.2f}\t{:10.2f}\t{:10.2f}\n".format(R_new[i],*pecs[i,:])
    txt4 += "{:.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\t{:10.2f}\
\t{:10.2f}\t{:10.2f}\n".format(R_new[i],*vv[i,:],*a_so[i,:])

with open("output_pecs.dat","w") as f:
    f.write(txt3)
with open("output_least_sq_fit.dat","w") as f:
    f.write(txt4)
