import numpy as np
import matplotlib.pyplot as plt
import glob

files = glob.glob("data/*.bin") # list of path+names for relevant files in directory
print(files)

plt.figure()
for f in files:
    plt.plot(np.fromfile(f).T[0,1], 'o')
plt.show()
