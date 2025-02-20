import numpy as np
import pandas as pd
from scipy.optimize import curve_fit as cf
import matplotlib.pyplot as plt;

df=pd.read_csv('C:\\Users\garet\Documents\MATLAB\GEPhysRepo\EL02\data2.csv')

f = df.freq.values
v2 = df.v2.values
v1 = df.v1.values

#print(df.head())

#print(f)

T2 = (v2 / v1)**2

def fit1(f,a,b,c):
    return 1/((1+c)**2+(b/a)**2*((f/a)**2+(a/f)**2-2))


dT2 = (2/(12**0.5))*(v2**4/v1**6+v2**2/v1**4)**0.5

p0 = 840, 3.6e4, 1.1#initial guesses
popt, pcov = cf(fit1, f, T2, p0,sigma=dT2)
xfit = np.arange(500,2000,1)
yfit = fit1(xfit, *popt)

res = T2 - fit1(f, *popt)
plt.subplot(211)
plt.scatter(f, T2, label='Raw Data')
plt.plot(xfit, yfit, 'r', label='Fitted Curve')
plt.legend()
plt.xscale('log')
plt.xlabel("Frequency / Hz")
plt.ylabel("Transmittance Squared")
plt.title("Transmittance Squared vs Frequency for LCR circuit with 1st Fit")

plt.annotate('a: ' + ('%.3e' % popt[0]) + 'b: ' + ('%.3e' % popt[1]) + 'c: ' + ('%.3e' % popt[2]) + '\nErrors (a, b, c): ' + ('%.3e' % np.sqrt(pcov[0][0]))+', '+ ('%.3e' % np.sqrt(pcov[1][1]))+', '+ ('%.3e' % np.sqrt(pcov[2][2])), [1000, 0.1])



plt.subplot(212)
plt.scatter(f, res, label = 'Residuals for 1st Fit')
plt.xscale('log')
plt.xlabel("Frequency / Hz")
plt.ylabel("Transmittance Error")
plt.title("1st Fit Residuals")

plt.subplots_adjust(hspace = 0.4)
plt.show()
