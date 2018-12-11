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

def graph(files, orderlist, name):
    plt.figure()
    for i in orderlist:
        plt.plot(
                np.fromfile(files[i])[::100], 'o-', label='L='+Llist[i], alpha=0.5
                   )

    plt.legend()
    plt.xlabel('Temperature')
    plt.ylabel(r"%s" % name)
    plt.grid()

Llist = ['40', '60', '80', '100']
#VarianceE= filelist("VarianceE")
#VarianceEorderlist = [1, 3, 2, 0]
#graph(VarianceE, VarianceEorderlist, "$C_V k_B / T^2 $")
#ExpectationEnergy = filelist("ExpectationEnergy")
#ExpectationEnergyOrderlist = [0, 1, 3, 2]
#graph(ExpectationEnergy, ExpectationEnergyOrderlist,"$\langle E \rangle$" )
#ExpectationMagnetism = filelist("ExpectationMagnetism")
#ExpectationMagnetismOrderlist = [1, 2, 3, 0]
#graph(ExpectationMagnetism, ExpectationMagnetismOrderlist, "$\langle |M| \rangle$")
#VarianceM = filelist("VarianceE")
#VarianceMorderlist = [1, 3, 2, 0]
#graph(VarianceM, VarianceMorderlist, "\chi k_B /T")

plt.show()
