import numpy as np
import matplotlib.pyplot as plt
import glob
import sys

assert( len(sys.argv) == 4 ), "Please provide a path, a temperature step and a run name to search for. Omit '.bin'"

path = str( sys.argv[1] )
TempStep = "Temp"+str( sys.argv[2] )
runName = str( sys.argv[3] )

files = glob.glob( path + "*"+ TempStep + runName + ".bin" ) # list of path+names for relevant files in directory
print(files)

cont = "y"
while cont == "y":
    plt.figure()
    filePos = int(input("what to plot?\n"))
    #cut1 = files[filePos].replace(path, "")
    #cut2 = cut1.replace(TempStep+runName+".bin", "")
    plottype = input("log plot?")
    plt.semilogx(np.fromfile(files[filePos])[1::], label=cut2)
    cont = input("continue(y/n)? ")
plt.grid()
plt.legend()
plt.xlabel('MC cycles')
plt.show()


runName =
variableName = input("name of variable to plot: ")

