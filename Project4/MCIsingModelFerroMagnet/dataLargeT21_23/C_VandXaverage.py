import numpy as np
import glob

C_vlist = glob.glob("VarianceE*")
X_list = glob.glob("VarianceM*")
print(C_vlist)
print(X_list)
Llist = ['40', '60', '80', '100']
Temp = np.fromfile("Temperature.bin")

CvTemp = np.zeros(len(Llist))
XTemp = np.zeros(len(Llist))


for i in range(len(Llist)):
    CvTemp[i] = Temp[np.argmax(np.fromfile(C_vlist[i]))]
    XTemp[i] = Temp[np.argmax(np.fromfile(X_list[i]))]


Tc_Cv = np.mean(CvTemp)
stderrCv = np.std(CvTemp)
Tc_X = np.mean(XTemp)
stderrX = np.std(XTemp)

print( "T_c from Cv: ", Tc_Cv, " +- ", stderrCv)
print( "T_c from X: ", Tc_X, " +- ", stderrX)

#C_vMean = np.zeros(len(C_vlist))
#C_vstderr = np.zeros(len(C_vlist))
#XMean = np.zeros(len(X_list))
#Xstderr = np.zeros(len(X_list))


#    C_vMean[i] = np.mean(Cv)
#    XMean[i] = np.mean(X)
#    C_vstderr[i] = np.std(Cv)
#    Xstderr[i] = np.std(X)
