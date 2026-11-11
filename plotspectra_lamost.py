import sys# JH 7/06/2024
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from PyAstronomy import pyasl
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Define the Gaussian function 
#def Gauss(x,X0, A, B, C): 
#    y = A*np.exp(-1*B*(x-X0)**2)+C 
#    a * (np.exp(-b*(x-x_0)**2))
#    return y 

df = pd.DataFrame()
df['fil1'] = None
df['flag'] = None

input1=sys.argv[1]
# LA LISTA DEBE TENER EXTENSION ".lis"
output1 = input1.replace('.lis', '_flag.csv')

with open(input1, 'r') as f:
    file1 = f.read().splitlines()
flag_ = []
fil1_ = []
for j in file1:
        j2 = j.split("/")
        with fits.open(j) as hdul:
          fl1 = hdul[2].data.field('FLUX')
          ivar1 = hdul[2].data.field('IVAR')
          wl2 = hdul[2].data.field('WAVELENGTH')
        fl = fl1[0,]
        efl=np.sqrt(1/ivar1[0,])
        wl1 = wl2[0,]
        wl = pyasl.vactoair2(wl1,mode='ciddor')

        plt.subplot(2,1,1)    
        plt.title(j2[len(j2) - 1])
        plt.plot(wl, fl)
        plt.fill_between(wl,fl-efl,fl+efl, alpha=0.5, color='red')
        
        plt.subplot(2,2,3)  ##Ha
        ind = np.where((wl>=6533) & (wl<=6593))
#        ind = np.where((wl>=7020) & (wl<=7100))# OIV
        wlHa = wl[ind]
        flHa = fl[ind]
        eflHa = efl[ind]
        plt.plot(wlHa, flHa)
        plt.fill_between(wlHa,flHa-eflHa,flHa+eflHa, alpha=0.3, color='red')
        plt.axvline(6563, linestyle = '--')

        plt.subplot(2,2,4)
#        ind = np.where((wl>=6700) & (wl<=6720))
        ind = np.where((wl>=6250) & (wl<=6350)) #  OI
        wlLi = wl[ind]
        flLi = fl[ind]
        eflLi = efl[ind]
#        parameters, covariance = curve_fit(Gauss, wlLi, flLi, [6707.0,20.0,1.5,750])
#        fit_X0 = parameters[0]
#        fit_A = parameters[1] 
#        fit_B = parameters[2]
#        fit_C = parameters[3]
#        fit_y = Gauss(wlLi, fit_X0, fit_A, fit_B, fit_C) 

        plt.plot(wlLi, flLi)
        plt.fill_between(wlLi,flLi-eflLi,flLi+eflLi, alpha=0.3, color='red')
#        plt.axvline(6705, linestyle = '--')
#        plt.axvline(6711, linestyle = '--')
#        plt.axvline(6698, linestyle = '--')
#        plt.axvline(6703, linestyle = '--')
        plt.axvline(6300, linestyle = '--')
#        plt.axvline(7068, linestyle = '--')
#        plt.plot(wlLi, fit_y, '-', label='fit')
        plt.show()
        print(" T Tauri: Y, N, or Q")
        T1=input()
        if T1 == "N" or T1 == "n" or T1=="Q" or T1=="q":
           print(j2[len(j2) - 1],T1)
        else:
           print(j2[len(j2) - 1],"Y")
           T1="Y"
        fil1_.append(j2[len(j2)-1])
        flag_.append(T1)
df['flag']=flag_
df['fil1']=fil1_
print(df)
df.to_csv(output1, index=False)
