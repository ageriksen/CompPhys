import matplotlib.pyplot as plt
import numpy as np

def read(path):
    collumn1 = []; collumn2 = []; collumn3 = []
    filename = open(path, 'r')
    filename.readline()
    filelines = filename.readlines()
    for line in filelines:
        words = line.split()
        collumn1.append(words[0])
        collumn2.append(words[1])
        collumn3.append(words[2])
        if (len(words) == 7):
            stepnumber = int(words[-1])
    print("number of steps: ", stepnumber)
    return stepnumber, np.array(collumn1), np.array(collumn2), np.array(collumn3)

    
def plot_bucklingbeam():
    return 0

def plot_1particleHO():
    return 0


def plot_2particleHO():
    return 0


if __name__ == '__main__':
    import sys
    step, state1, state2, state3 = read(sys.argv[1])
