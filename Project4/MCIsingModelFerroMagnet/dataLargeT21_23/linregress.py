import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy import stats


def fit(x, y):
    # linear regression of arrays x and y
    slope = stats.linregress(x,y)[0]
    intercept = stats.linregress(x,y)[1]
    std_err = stats.linregress(x, y)[4]
    return slope,intercept, std_err

def filelist(Lval):
    #picks out a list of files with L values equivalent to Lval
    files = glob.glob("VarianceE*L"+Lval+".bin")
    print( files )
    return files

def maximum(array):
    #Extracting the indexes of the maximal values in "array"
    return np.argmax(array)

C040 = np.fromfile(filelist('040')[0])
C060 = np.fromfile(filelist('060')[0])
C080 = np.fromfile(filelist('080')[0])
C100 = np.fromfile(filelist('100')[0])
Temp = np.fromfile("Temperature.bin")
max040 = maximum(C040)
max060 = maximum(C060)
max080 = maximum(C080)
max100 = maximum(C100)
maxarray = np.array((
    C040[max040], C060[max060], C080[max080], C100[max100]
    ))
Temps = np.array((
    Temp[max040],Temp[max060],Temp[max080],Temp[max100]
    ))
print( maxarray )
print( Temps )
slope, intercept, std_err = fit(maxarray, Temps )
print( "==========================================")
print( " T_C( L ) = T_C( L=inf ) + a/L \n For L=inf, a/L = 0 \n T( L=inf ) approx: ", intercept )
print( "rel_error: ", std_err)
L = np.arange(40, 101, 20)
print( L )
plt.plot(L, L*slope + intercept, label=r"$ T_C ( L ) $" )
plt.plot(L, Temps, 'o', label=r"$T(C_{V,max} )$")
plt.xlabel(r"$L$")
plt.ylabel(r"$T_C(L)$")
plt.grid()
plt.legend()
plt.show()
