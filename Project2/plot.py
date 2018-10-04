import matplotlib.pyplot as plt
import numpy as np

def read2particle(path):
    collumn1 = []
    filename = open(path, 'r')
    filename.readline()
    filelines = filename.readlines()
    for line in filelines:
        words = line.split()
        collumn1.append(float(words[0]))
        if (len(words) == 7):
            stepnumber = int(words[-1])
    print("number of steps: ", stepnumber)
    return stepnumber, np.array(collumn1)

    
def plot_bucklingbeam():
    return 0

def plot_1particleHO():
    return 0


def plot_2particleHO():
    number = np.zeros(4)
    path = []
    for i in range(4):
        print( i )
        path.append( "omega"+str(i)+"dim200Eigval.txt")
    number[0], stateomega0 = read2particle(path[0])
    number[1], stateomega1 = read2particle(path[1])
    number[2], stateomega2 = read2particle(path[2])
    number[3], stateomega3 = read2particle(path[3])
    step = 10./number
    prob1 = stateomega0**2; prob2 = stateomega1**2; prob3 = stateomega2**2; prob4 = stateomega3**2;
    rho = np.linspace(0, 10, number[0])
    plt.figure()
    plt.plot(rho, prob1)
    plt.plot(rho, prob2)
    plt.plot(rho, prob3)
    plt.plot(rho, prob4)
    plt.xlabel(r"$\rho$")
    plt.ylabel(r"$|\psi |^2$")
    plt.show()

    return 0


if __name__ == '__main__':
    import sys
    plot_2particleHO()
