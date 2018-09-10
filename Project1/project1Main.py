"""
This is a program to execute and link project 1 in computational pysics. 
the main task is to solve a simple differential equation
u¨ = F(x) | u[i+1] + u[i-1] - 2u[i] = h*f(x[i]) 
with discretization and approximation through a Taylor expansion. 
"""
import sys
import numpy as np
from numba import jit
import matplotlib.pyplot as plt
#importing project modules
import TridiagAlgo as tri

#declare variables
n = np.zeros(int(sys.argv[1]), dtype=np.int32)
h = np.zeros(len(n))

#filling arrays
for i in range(1, len(n)+1):
    n[i-1] = 1*10**i
    h[i-1] = float(1)/(n[i-1]+1)

#positional- and function-arrays
for i in range(len(n)):
    print('length of arrays ~ ', n[i])
    x = np.linspace(0,1,n[i]+2)
    u = np.zeros(int(n[i]+2))
    #exact solution
    f = 100*np.exp(-10*x)
    d = (h[i]**2)*f
    #constructing tridiagonal matrix A
    a = -1*np.ones(int(n[i]))
    b = 2*np.ones(int(n[i]+1))
    c = -1*np.ones(int(n[i]))
    # calling tri module and solving for u
    dcomp = tri.TriSubstitution(u)
    dcomp(a, b, c, d)
    u, FLOPS = dcomp.Substitute(n[i])
    
    u_ex = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)
    
    #ploting results
    plt.figure()
    plt.plot(x, u, label='numerical')
    plt.plot(x, u_ex, label='exact')
    plt.legend()
    plt.title('numerical and analytical solution of u''=f(x), n = '+str(n[i]))
    plt.xlabel('position x')
    plt.ylabel('value solution u')
    plt.show()
