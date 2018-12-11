import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

variables = [
    r"$\langle |M| \rangle / N_{spins}$",
    r"$ \chi / N_{spins}$",
    r"$\langle E \rangle / N_{spins}$",
    r"$ C_V / N_{spins}$"
    ]
temp = ['10', '24']

files1 = glob.glob("Expectation*RandomL20T"+temp[0]+".bin")
files2 = glob.glob("Expectation*RandomL20T"+temp[1]+".bin")
files3 = glob.glob("Expectation*UpL20T"+temp[0]+".bin")
files4 = glob.glob("Expectation*UpL20T"+temp[1]+".bin")
print(files1)
print(files2)
print(files3)
print(files4)

def graph(f, d, g):
    plt.figure()
    plt.semilogx(np.fromfile(files1[f])[::1000], 'o-', label='T='+temp[0]+" random", alpha=0.4)
    plt.semilogx(np.fromfile(files2[d])[::1000], 'o-', label='T='+temp[1]+" random", alpha=0.4)
    plt.semilogx(np.fromfile(files3[f])[::1000], 'o-', label='T='+temp[0]+" ordered", alpha=0.4)
    plt.semilogx(np.fromfile(files4[f])[::1000], 'o-', label='T='+temp[1]+" ordered", alpha=0.4)
    plt.legend()
    plt.ylabel('scaled with coupling J=1')
    plt.xlabel('MC cycles')
    plt.grid()

graph(0, 1, 0)
graph(1, 0, 2)
plt.show()
