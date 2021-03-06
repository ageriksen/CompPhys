import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def readfile(path, dim):
    print("read file: ",path)
    if dim == 2:
        return np.loadtxt(path, usecols=(0,1), unpack=True)
    elif dim == 3:
        return np.loadtxt(path, unpack=True)

def plotF3D(nr, paths, graphlabels, plottitle, quarter=False):
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_zlabel(r"$z$")  
    for i in range(nr):
        ax.plot(*readfile(path[i], 3), label=graphlabels[i])
    ax.legend()
    #plt.title(plottitle)
    plt.show()
    return 0
    
def plotF2D(nr, paths, graphlabels, plottitle, quarter=False):
    fig = plt.figure()
    plt.xlabel(r"$x$")
    plt.ylabel(r"$y$")
    for i in range(nr):
        plt.plot(*readfile(path[i], 2), label=graphlabels[i])
    plt.legend()
    plt.title(plottitle)
    plt.show()
    return 0


def subplotF3D(nr, path1, path2, path3, path4, labels):
    fig = plt.figure()

    ax = fig.add_subplot(2, 2, 1, projection='3d')
    #ax = fig.gca(projection="3d")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_zlabel(r"$z$")  
    for i in range(nr):
        ax.plot(*readfile(path1[i], 3), label=labels[i])
    ax.legend()
    
    ax = fig.add_subplot(2, 2, 2, projection='3d')
    #ax = fig.gca(projection="3d")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_zlabel(r"$z$")  
    for i in range(nr):
        ax.plot(*readfile(path2[i], 3), label=labels[i])
    ax.legend()
   
    ax = fig.add_subplot(2, 2, 3, projection='3d')
    #ax = fig.gca(projection="3d")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_zlabel(r"$z$")  
    for i in range(nr):
        x, y, z = readfile(path3[i], 3)
        ax.plot(x[:(int(len(x)/6))], y[:(int(len(y)/6))], z[:(int(len(z)/6))], label=labels[i])
    ax.legend()
    
    ax = fig.add_subplot(2, 2, 4, projection='3d')
    #ax = fig.gca(projection="3d")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$y$")
    ax.set_zlabel(r"$z$")  
    for i in range(nr):
        x, y, z = readfile(path3[i], 3)
        ax.plot(x[:(int(len(x)/6))], y[:(int(len(y)/6))], z[:(int(len(z)/6))], label=labels[i])
    ax.legend()
    #plt.title(plottitle)
    plt.show()
    return 0

if __name__ == "__main__":
    import sys
#    path = ["mobileSun3Bodysun.dat", "mobileSun3Bodyearth.dat", "mobileSun3Bodyjupiter.dat"] 
#    labels = ["sun", "earth", "jupiter"]
#    path = [
#            "fullSystemsun.dat",
#            "fullSystemmercury.dat",
#            "fullSystemvenus.dat",
#            "fullSystemearth.dat",
#            "fullSystemmars.dat",
#            "fullSystemjupiter.dat",
#            "fullSystemsaturn.dat",
#            "fullSystemuranus.dat",
#            "fullSystemneptune.dat"
#            ]
#    labels = [
#            "Sun",
#            "Mercury",
#            "Venus",
#            "Earth",
#            "Mars",
#            "Jupiter",
#            "Saturn",
#            "Uranus",
#            "Neptune"
#            ]
#    path = ["relativisticPrecessionmercury.dat"]
#    labels = ["Mercury"]
#    nr = int(len(path)) # Department of redundancy department
#    title=r"Sun Earth Jupiter system"
#    plotF2D(nr, path, labels, title)
#    plotF3D(nr, path, labels, title)
    nr = 2
    path1 = ["stationarySun3Bodyearth.dat", 
            "stationarySun3Bodyjupiter.dat"]
    path2 = ["10xJupiter3Bodyearth.dat", 
            "10xJupiter3Bodyjupiter.dat"]
    path3 = ["1kxJupiter3Bodyearth.dat",
            "1kxJupiter3Bodyjupiter.dat"]
    path4 = ["1kxJupiter3Bodyearth.dat",
            "1kxJupiter3Bodyjupiter.dat"]
    labels = ["earth", "jupiter"]
    subplotF3D(nr, path1, path2, path3, path4, labels)
