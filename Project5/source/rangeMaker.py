"""
short script to make ranges of alpha and beta
"""
import numpy as np

def makeAlpha(start, stop, step):
    alphaRange = np.arange(start, stop, step)
    np.savetxt("resources/alpha.dat", alphaRange, fmt='%.2f')

def makeBeta(start, stop, step):
    betaRange = np.arange(start, stop, step)
    np.savetxt("resources/beta.dat", betaRange, fmt='%.2f')

# initial arays:
makeAlpha(0.8, 1.23, 0.05)
makeBeta(0.8, 1.23, 0.05)
