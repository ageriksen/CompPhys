import numpy as np
import matplotlib.pyplot as plt
import sys
import glob

"""
SO Much redundancy here, but I need to just crank it out for now.
"""

def getData( fileName ):
    return np.loadtxt(fileName, skiprows=1)

def getLabels( fileName ):
    f = open(fileName)
    line = f.readline()
    return line.split()

def getValues( dataMatrix, i ):
    return dataMatrix[:, i]

def minColumn( dataMatrix, i ):
    return np.argmin( dataMatrix[:,i] )

def cutElements( dataMatrix, controlColumn ):
    indexList = []
    column = dataMatrix[:, controlColumn]
    print(column)
    previous = column[0]
    for element in column :
        if( element != previous ):
            print(np.where(column==element) )
            indexList.append(np.where(column==element))
        prev = element
    return indexList


# plotting energy vs. Alphavalues
fileName = "../data/trial1Full/trialwf1fullOmega1.000000.dat"
labels = getLabels( fileName )
data = getData( fileName )
#cutIndices = cutElements( data, 0 )
#print(cutIndices)
alpha = getValues( data, 0 )
energy = getValues( data, 1 )
variance = getValues( data, 2 )

plt.figure()
plt.errorbar(alpha, energy, yerr=variance )
#plt.plot( getValues( datam 1 ), getValues( data, 2 )
plt.xlabel(labels[0])
plt.ylabel(labels[1])
plt.show()
