import numpy as np
import matplotlib.pyplot as plt
import glob
import sys


files = glob.glob( input("which folder?"))
print(files)

#plt.figure()
plt.semilogx(
    np.fromfile(files[input("which index to plot?")] ),
    alpha=0.7,
    label=r""+input("label?")
    )
plt.grid()
plt.legend()
plt.show()
