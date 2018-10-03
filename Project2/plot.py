import matplotlib.pyplot as plt
import numpy as np

def read2particle(path):
    collumn1 = []
    filename = open(path, 'r')
    filename.readline()
    filelines = filename.readlines()
    for line in filelines:
        words = line.split()
        collumn1.append(words[0])
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
    number[1], stateomega0 = read2particle(path[1])
    number[2], stateomega0 = read2particle(path[2])
    number[3], stateomega0 = read2particle(path[3])
    step = 10./number
    

    return 0


if __name__ == '__main__':
    import sys
    plot_2particleHO()
