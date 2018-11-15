import numpy as np
import matplotlib.pyplot as plt
import glob

files = glob.glob("data/*.bin") # list of path+names for relevant files in directory

eV = np.zeros((len(files),2))
for i in range(len(files)):
    expVal = np.fromfile(files[i])
    eV[i,0], eV[i,1] = expVal[0], expVal[1]

print(eV)

plt.figure()
plt.plot(eV[:,1], eV[:,0], 'o')
plt.show()
