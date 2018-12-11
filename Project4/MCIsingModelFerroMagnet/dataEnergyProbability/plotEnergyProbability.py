import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import glob
import sys

variables = [
    r"$\langle |M| \rangle / N_{spins}$",
    r"$ \chi / N_{spins}$",
    r"$\langle E \rangle / N_{spins}$",
    r"$ C_V / N_{spins}$"
    ]
temp = ['10', '24']

files = glob.glob("Energies*")
print(files)
Energies10 = np.fromfile(files[1])[::1000]
Energies24 = np.fromfile(files[0])[::1000]

plt.subplot(211)
plt.hist(Energies10, density=True)
plt.ylabel(r"$\Pr{(E)}$")
plt.subplot(212)
plt.hist(Energies24, density=True)
plt.xlabel('Energies')
plt.ylabel(r"$\Pr{(E)}$")
plt.show()
