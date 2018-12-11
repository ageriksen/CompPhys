import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

def filelist(string):
    files = glob.glob(string+"*")
    print(files)
    return files

def graph(files, orderlist, name, T, pos):
    plt.subplot(pos)
    plt.semilogx(
            files[::10000],
            'o-',
            label='L='+T, alpha=0.5
               )

    plt.legend()
    plt.xlabel(r'$T\cdot 10^2$')
    plt.ylabel(r"%s" % name)
    plt.grid()

Tlist = ['1.0', '2.4']
AcceptedList = filelist("Accepted*")
AcceptedOrderlist = [0, 1]
Accepted1 = np.fromfile(AcceptedList[0])
Accepted2 = np.fromfile(AcceptedList[1])
Accepted1 = Accepted1/np.arange(len(Accepted1))
Accepted2 = Accepted2/np.arange(len(Accepted2))

graph(Accepted1, 0, "", '1.0̈́', 211)
graph(Accepted2, 1, "accepted changes", '2.4', 212)

plt.show()
