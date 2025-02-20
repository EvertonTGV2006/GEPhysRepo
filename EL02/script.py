import numpy as np
import pandas as pd
from scipy.optimize import curve_fit as cf
import matplotlib.pyplot as plt;

df=pd.read_csv('C:\\Users\garet\Documents\MATLAB\GEPhysRepo\EL02\data.csv')

tf = df.freq.values
tv2 = df.v2.values
tv1 = df.v1.values

f = np.resize(tf, tf.size - 2)
v2 = np.resize(tv2, tv2.size - 2)
v1 = np.resize(tv1, tv1.size - 2)
#print(f)

T2 = (v2 / v1)**2

def fit1(f,a,c):
    return 1/(a+c*f**2)


dT2 = (2/(12**0.5))*(v2**4/v1**6+v2**2/v1**4)**0.5

p0 = 1, 1.3e-6 #initial guesses
popt, pcov = cf(fit1, f, T2, p0,sigma=dT2)

xfit = np.arange(100,10000,1)
yfit = fit1(xfit, popt[0], popt[1])

res = T2 - fit1(f, popt[0], popt[1])
plt.subplot(321)
plt.scatter(f, T2, label='Raw Data')
plt.plot(xfit, yfit, 'r', label='Fitted Curve')
plt.legend()
plt.xscale('log')
plt.xlabel("Frequency / Hz")
plt.ylabel("Transmittance Squared")
plt.title("Transmittance Squared vs Frequency for RC circuit with 1st Fit")

plt.annotate('a: ' + ('%.3e' % popt[0]) + 'c: ' + ('%.3e' % popt[1]) + '\nErrors (a, c): ' + ('%.3e' % np.sqrt(pcov[0][0]))+', '+ ('%.3e' % np.sqrt(pcov[1][1])), [100, 0])



plt.subplot(322)
plt.scatter(f, res, label = 'Residuals for 1st Fit')
plt.xscale('log')
plt.xlabel("Frequency / Hz")
plt.ylabel("Transmittance Error")
plt.title("1st Fit Residuals")


def fit2(f, a, b, c):
    return 1/(a + b*f + c*f**2)

p1 = 1, 1e-6, 1e-6
p2opt, p2cov = cf(fit2, f, T2, p1, sigma=dT2)

yfit2 = fit2(xfit, p2opt[0], p2opt[1], p2opt[2])
res2 = T2 - fit2(f, *p2opt)

plt.subplot(323)
plt.scatter(f, T2, label = 'Raw Data')
plt.plot(xfit, yfit2,'r', label = '2nd Improved Fit')
plt.legend()
plt.xscale('log')
plt.xlabel('Frequency / Hz')
plt.ylabel('Transmittance Squared')
plt.title('2nd Fit')

plt.annotate('a: ' + ('%.3e' % p2opt[0]) + 'b: ' + ('%.3e' % p2opt[1]) + 'c: ' + ('%.3e' % p2opt[2]) + '\nErrors (a, b, c): ' + ('%.3e' % np.sqrt(p2cov[0][0]))+', '+ ('%.3e' % np.sqrt(p2cov[1][1]))+', '+ ('%.3e' % np.sqrt(p2cov[2][2])), [100, 0])

plt.subplot(324)
plt.scatter(f, res2)
plt.xscale('log')
plt.xlabel('Frequency / Hz')
plt.ylabel('Transmittance Error')
plt.title('Residuals on the 2nd Fit')

def fit3(f, a, b, c, d, e):
    return (a + b*f + c*f**2 + d*f**3 + e*f**4)**-1

p2 = 1, 1e-6, 1e-6, 1e-6, 1e-6
p3opt, p3cov = cf(fit3, f, T2, p2, sigma=dT2)

yfit3 = fit3(xfit, *p3opt)
res3 = T2 - fit3(f, *p3opt)

plt.subplot(325)
plt.scatter(f, T2, label = 'Raw Data')
plt.plot(xfit, yfit3,'r', label = '2nd Improved Fit')
plt.legend()
plt.xscale('log')
plt.xlabel('Frequency / Hz')
plt.ylabel('Transmittance Squared')
plt.title('3rd Fit')

plt.annotate('a: ' + ('%.3e' % p3opt[0]) + 'b: ' + ('%.3e' % p3opt[1]) + 'c: ' + ('%.3e' % p3opt[2]) + 'd: ' + ('%.3e' % p3opt[3]) +'e: ' + ('%.3e' % p3opt[4]) + '\nErrors (a, b, c, d, e): ' + ('%.3e' % np.sqrt(p3cov[0][0]))+', '+ ('%.3e' % np.sqrt(p3cov[1][1]))+', '+ ('%.3e' % np.sqrt(p3cov[2][2])) + ', '+ ('%.3e' % np.sqrt(p3cov[3][3])) + ', '+ ('%.3e' % np.sqrt(p3cov[4][4])), [100, 0])

plt.subplot(326)
plt.scatter(f, res3)
plt.xscale('log')
plt.xlabel('Frequency / Hz')
plt.ylabel('Transmittance Error')
plt.title('Residuals on the 3rd Fit')

plt.subplots_adjust(top = 0.95, bottom=0.1, hspace = 0.5)
plt.show()
