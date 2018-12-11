"""
Plot <E>, <|M|>, Cv and X as func of T
calibrated by first plotting the file lists, to
get the correct indices, as glob seems to pick
files according to memory order rather than alphabetical.
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import sys


def filelist(string):
    files = glob.glob(string+"*")
    print(files)
    return files

def graph(files, orderlist, name, b):
    plt.subplot(b)
    for i in orderlist:
        f = files[i]
        print(len(f))
        plt.plot(
                T,
                np.fromfile( f),
                'o-',
                label='L='+Llist[i],
                alpha=0.5
                   )

    plt.legend()
    plt.xlabel('Temperature')
    plt.ylabel(r"%s" % name)
    plt.grid()


#=====================

plt.figure()
TempList = filelist("Temp")[0]
T = np.fromfile(TempList)
Llist = ['40', '60', '80', '100']

ExpectationEnergy = filelist("ExpectationEnergy")
ExpectationEnergyOrderlist = [0, 1, 3, 2]
graph(ExpectationEnergy, ExpectationEnergyOrderlist, r"$\langle E \rangle$", 221)

VarianceE= filelist("VarianceE")
VarianceEorderlist = [1, 3, 2, 0]
graph(VarianceE, VarianceEorderlist, r"$C_V/k_B t^2$", 223)

ExpectationMagnetism = filelist("ExpectationMagnetism")
ExpectationMagnetismOrderlist = [1, 2, 3, 0]
graph(ExpectationMagnetism, ExpectationMagnetismOrderlist, r"$\langle |M| \rangle $", 222)

VarianceM = filelist("VarianceM")
VarianceMorderlist = [0, 1, 3, 2]
graph(VarianceM, VarianceMorderlist, r"$\chi/k_B T$", 224)

plt.show()
