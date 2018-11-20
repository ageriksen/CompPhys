import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

path = input()
LandT = input()

E_files = glob.glob( path+"ExpectationEnergy"+LandT )
M_files = glob.glob( path+"ExpectationMagnetism"+LandT )
VarE_files = glob.glob( path+"VarianceE"+LandT )
VarM_files = glob.glob( path+"VarianceM"+LandT )
print(E_files)
print(M_files)
print(VarE_files)
print(VarM_files)

plt.figure()
plt.semilogx(
   np.fromfile(files[E_files])[2::],
   alpha=0.7,
   label=r"\langle E \rangle / N"
   )
plt.grid()

plt.figure()
plt.semilogx(
       np.fromfile(files[M_files])[2::],
       alpha=0.7,
       label=r"\langle |M| \rangle / N"
       )
plt.grid()

plt.figure()
plt.semilogx(
       np.fromfile(files[VarE_files])[2::],
       alpha=0.7,
       label=r"\sigma_E \def C_V\cdot k_B"
       )
plt.grid()

plt.figure()
plt.semilogx(
       np.fromfile(files[VarM_files])[2::],
       alpha=0.7,
       label=r"\sigma_{|M|} \def \chi_{|M|}\cdot k_B"
       )
plt.grid()

def Energy():
    for i in range(len(L_list)):
        for j in range(len(T_list)):

#plt.show()
